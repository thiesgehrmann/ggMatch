
import random
import inspect, os
__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
__PC_DIR__ = "%s/pipeline_components" % __INSTALL_DIR__

###############################################################################

import json
dconfig = json.load(open("%s/defaults.json"% __PC_DIR__, "r"))
dconfig.update(config)

###############################################################################

from pipeline_components import blastutils as butils
from pipeline_components import utils as utils



###############################################################################

__RUN_DIR__ = os.path.abspath(dconfig["outdir"]) + "/run"
__BLASTDB_OUTDIR__ = "%s/blastdb" % __RUN_DIR__
__ITERATION_OUTDIR__ = "%s/iterations" % __RUN_DIR__
__FINAL_OUTDIR__ = "%s/final"% __RUN_DIR__
__GRAPH_OUTDIR__ = "%s/graph" % __RUN_DIR__
__CLUSTALO_OUTDIR__ = "%s/alignment" % __RUN_DIR__
__COMPARE_OUTPUT__ = "%s/compare" % __RUN_DIR__
__PHYLO_OUTPUT__ = "%s/phylogeny" % __RUN_DIR__
__RECIPROCALBLAST_OUTPUT__ = "%s/reciprocalBlasts" % __RUN_DIR__
__VALIDATION_OUTDIR__ = "%s/validation" % __RUN_DIR__

###############################################################################

querySetProts = "%s/querySet.fasta" % __ITERATION_OUTDIR__

dconfig["genomes"][dconfig["querySetName"]] = { "prots" : querySetProts }

###############################################################################

rule makeBlastDB:
  input:
    proteins = lambda wildcards: dconfig["genomes"][wildcards.genome]["prots"]
  output:
    db = "%s/{genome}.blastdb" % __BLASTDB_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    #makeblastdb -dbtype prot -in {input.proteins} -out {output.db}
    diamond makedb --in {input.proteins} -d {output.db}
    touch {output.db}
  """

rule blastDBs:
  input:
    dbs = expand("%s/{genome}.blastdb" % __BLASTDB_OUTDIR__, genome=sorted(dconfig["genomes"].keys()))

###############################################################################

rule genIterationJSON:
  input:
    queries = lambda wildcards: "%s/%s/iteration_%d.fasta" % (__ITERATION_OUTDIR__, wildcards.query, int(wildcards.iteration)-1)
  output:
    json = "%s/{query}/{iteration}.json" % (__ITERATION_OUTDIR__)
  run:
    query = wildcards.query
    iteration = int(wildcards.iteration)

    speciesPresent = [ x.split(":")[0] for x in utils.loadFasta(input.queries) ]

    blastDBs = { genome: "%s/%s.blastdb" % (__BLASTDB_OUTDIR__, genome) for genome in sorted(dconfig["genomes"].keys()) if genome not in speciesPresent  }
    prots    = { genome: config["genomes"][genome]["prots"] for genome in sorted(dconfig["genomes"].keys()) if genome not in speciesPresent}

    needNext = True
    if len(prots) == 0:
      needNext = False
    #fi

    if needNext and (int(wildcards.iteration) > 1):
      prevPrevQueryFile = "%s/%s/iteration_%d.fasta" % (__ITERATION_OUTDIR__, wildcards.query, int(wildcards.iteration)-2)
      prevPrevSpeciesPresent = [ x.split(":")[0] for x in utils.loadFasta(prevPrevQueryFile) ]

      if set(speciesPresent) == set(prevPrevSpeciesPresent):
        needNext = False
      #fi
    #fi

    C = { "blastdbs" : { genome: "%s/%s.blastdb" % (__BLASTDB_OUTDIR__, genome) for genome in sorted(dconfig["genomes"].keys()) if genome not in speciesPresent  },
          "prots" : { genome: config["genomes"][genome]["prots"] for genome in sorted(dconfig["genomes"].keys()) if genome not in speciesPresent},
          "allblastdbs" : { genome: "%s/%s.blastdb" % (__BLASTDB_OUTDIR__, genome) for genome in sorted(dconfig["genomes"].keys())  },
          "allprots" : { genome: config["genomes"][genome]["prots"] for genome in sorted(dconfig["genomes"].keys()) },
          "queryName" : wildcards.query,
          "iteration" : int(wildcards.iteration),
          "queries" : input.queries,
          "reciprocalDir" : __RECIPROCALBLAST_OUTPUT__,
          "outdir" : "%s/%s/iteration_%d.run" % ( __ITERATION_OUTDIR__, query, iteration),
          "final" : "%s/%s/iteration_%d.fasta" % ( __ITERATION_OUTDIR__, query, iteration),
          "matchedList" : "%s/allMatched.tsv" % (__ITERATION_OUTDIR__),
          "querySetName" : dconfig["querySetName"],
          "evalue" : dconfig["evalue"],
          "needNext" : needNext,
          "pcDir" : __PC_DIR__ }
    with open(output.json, "w") as ofd:
      ofd.write(json.dumps(C, indent=4))

rule initIteration0:
  input:
    query = lambda wildcards: dconfig["queries"][wildcards.query]
  output:
    query = "%s/{query}/iteration_0.fasta" % (__ITERATION_OUTDIR__)
  run:
    F  = utils.loadFasta(input.query)
    NF = {}
    for ident in F:
      genome = ident.split(':')[0]
      if genome not in dconfig["genomes"]:
        NF["%s:%s" % (dconfig["querySetName"], ident)] = F[ident]
      else:
        NF[ident] = F[ident]
      #fi
    #efor
    utils.writeFasta(NF.items(), output.query)

    with open("%s/allMatched.tsv" % __ITERATION_OUTDIR__, "a+") as ofd:
      for seed in NF:
        ofd.write("%s\t%d\tQUERY:%s\t%s\t%d\n" % (wildcards.query, 0, wildcards.query, seed, 1))

rule genQuerySetProts:
  input:
    files = [ dconfig["queries"][query] for query in dconfig["queries"].keys() ]
  output:
    querySetProts = dconfig["genomes"][dconfig["querySetName"]]["prots"]
  run:
    F = {}
    for f in input.files:
      F.update(utils.loadFasta(f))
    #efor
    utils.writeFasta(F.items(), output.querySetProts)
      

rule runIteration:
  input:
    blastdbs = rules.blastDBs.input.dbs,
    prevQuery = lambda wildcards: "%s/%s/iteration_%s.fasta" % (__ITERATION_OUTDIR__, wildcards.query, int(wildcards.iteration)-1),
    json = lambda wildcards: "%s/%s/%s.json" % (__ITERATION_OUTDIR__, wildcards.query, int(wildcards.iteration))
  output:
    nextQuery = "%s/{query}/iteration_{iteration,[1-9]}.fasta" % __ITERATION_OUTDIR__
  threads: 99 # Make sure that only one of these is run at a time!
  params:
    pc_dir = __PC_DIR__,
    matchedList = "%s/allMatched.tsv" % __ITERATION_OUTDIR__
  shell: """
    echo "####################################################################"
    echo "# Query: {wildcards.query}"
    echo "# Iteration: {wildcards.iteration}"
    echo "####################################################################"

    touch "{params.matchedList}"
    need_next=`cat "{input.json}" | grep "needNext" | cut -d ':' -f2 | tr -d ' ,'`
    if [ "$need_next" == "false" ]; then
      ln -s "{input.prevQuery}" "{output.nextQuery}"
    else
      {params.pc_dir}/run_iteration.sh {input.json} {threads} > /dev/null
    fi
  """
    
#rule runAllIterations:
#  input:
#    groups = lambda wildcards: expand("%s/%s/iteration_{iteration}.fasta" % (__ITERATION_OUTDIR__, wildcards.query), iteration=[ str(x) for x in range(1,10)])

rule moveFinal:
  input:
    groups = lambda wildcards: expand("%s/%s/iteration_{iteration}.fasta" % (__ITERATION_OUTDIR__, wildcards.query), iteration=[ str(x) for x in range(1,dconfig["max_iterations"])]),
    lastIter = lambda wildcards: "%s/%s/iteration_%d.fasta" % (__ITERATION_OUTDIR__, wildcards.query, dconfig["max_iterations"]-1),
  output:
    finalLoc = "%s/{query}.fasta" % (__FINAL_OUTDIR__),
  run:
    F = [ (ident, seq) for (ident, seq) in utils.loadFasta(input.lastIter).items() if "__INPUTQUERY__" not in ident ]
    utils.writeFasta(F, output.finalLoc)

rule final:
  input:
    final = expand("%s/{query}.fasta" % __FINAL_OUTDIR__, query=dconfig["queries"].keys()),

###############################################################################

rule matchedListToGraph:
  input:
    files = rules.final.input.final
  output:
    nodesCSV = "%s/nodes.csv" % __GRAPH_OUTDIR__,
    edgesCSV = "%s/edges.csv" % __GRAPH_OUTDIR__
  run:
    ml = utils.readColumnFile("%s/allMatched.tsv" % __ITERATION_OUTDIR__, "query iteration origin target keep", types="str int str str int")

    nodes = {}
    nodeAttr = {}
    edges = {}
    nodeIDcounter = 1
    for match in ml[::-1]:
      if match.keep == 0:
        continue
      #fi

      originGenome, originGene = utils.splitgg(match.origin)
      targetGenome, targetGene = utils.splitgg(match.target)

      origin = match.origin
      target = match.target
      query = match.query
      iteration = match.iteration

      if query not in nodes:
        nodes[query] = {}
        nodeAttr[query] = {}
        edges[query] = []
      #fi
      if origin not in nodes[query]:
        nodes[query][origin] = "node%d" % nodeIDcounter
        nodeIDcounter += 1
        nodeAttr[query][origin] = { "iteration": iteration-1, "label": "\"%s\"" % origin} #=\"%s\"]" % (iteration, origin)
      #fi

      if target not in nodes[query]:
        nodes[query][target] = "node%d" % nodeIDcounter
        nodeIDcounter += 1
        nodeAttr[query][target] = {"iteration": iteration, "label": "\"%s\"" % target} #"[peripheries=%d,label=\"%s\"]" % (iteration, target)
      #fi

      edges[query].append( (origin, target, match.keep))

    #efor


    # Write graph data in csv format, for gephi
    with open(output.nodesCSV, "w") as ofd:
      ofd.write("Id,query,label,iteration\n")
      for query in nodes:
        for node in nodes[query]:
          ofd.write("%s,%s,%s,%s\n" % (nodes[query][node], query, nodeAttr[query][node]["label"], nodeAttr[query][node]["iteration"]))
        #efor
      #efor
    #ewith

    with open(output.edgesCSV, "w") as ofd:
      ofd.write("source,target,validated\n")
      for query in edges:
        for (origin, target, keep) in edges[query]:
          ofd.write("%s,%s,%d\n" % (nodes[query][origin], nodes[query][target], keep))
        #efor
      #efor
    #ewith

###############################################################################

rule alignGroup:
  input:
    seqs = lambda wildcards: "%s/%s.fasta" % (__FINAL_OUTDIR__, wildcards.query)
  output:
    aln = "%s/{query}.aln.msf" % __CLUSTALO_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    clustalo -i "{input.seqs}" --threads 1 -o "{output.aln}"
  """

rule alignedPrettyPlot:
  input:
    aln = lambda wildcards: "%s/%s.aln" % (__CLUSTALO_OUTDIR__, wildcards.query)
  output:
    plot = "%s/{query}.pdf" % __CLUSTALO_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    prettyplot 
  """

###############################################################################

rule alignToInitial:
  input:
    initial = lambda wildcards: dconfig["queries"][wildcards.query],
    final   = lambda wildcards: "%s/%s.fasta" % (__FINAL_OUTDIR__, wildcards.query)
  output:
    res = "%s/blasts/cmp.{query}.tsv" % __COMPARE_OUTPUT__
  conda: "%s/conda.yaml" % __PC_DIR__
  threads: 4
  params:
    blastfields = butils.blastfields_positive
  shell: """
    makeblastdb -dbtype prot -in {input.initial} -out {output.res}.db
    blastp -query {input.final} -db {output.res}.db -outfmt "6 {params.blastfields}" -out {output.res} -num_threads {threads}
  """

rule compareToInitial:
  input:
    bres = expand("%s/blasts/cmp.{query}.tsv" % __COMPARE_OUTPUT__, query=sorted(dconfig["queries"].keys()))
  output:
    cmpTable = "%s/cmpTable.tsv" % __COMPARE_OUTPUT__
  run:
    cmps = {}
    genomes = [ g for g in dconfig["genomes"] if g != dconfig["querySetName"] ]
    for query, bresFile in zip(sorted(dconfig["queries"].keys()), input.bres):
      cmps[query] = {}
      bres = utils.indexListBy(butils.readBlastFile(bresFile, butils.blastfields_positive), lambda x: utils.splitgg(x.qseqid)[0])
      for genome in genomes:
        if genome not in bres:
          cmps[query][genome] = -1
        else:
          bestHit = butils.bestHit(bres[genome])
          score = float(bestHit.positive) / float(bestHit.qlen)
          cmps[query][genome] = int(score * 100)
        #fi
      #efor
    #efor

    with open(output.cmpTable, "w") as ofd:
      ofd.write("Genome\t%s\n" % "\t".join(sorted(dconfig["queries"].keys())))
      for genome in genomes:
        ofd.write("%s\t%s\n" % (genome, "\t".join( "%d" % cmps[query][genome] for query in sorted(dconfig["queries"].keys()))))
      #efor
    #ewith
          

###############################################################################
###############################################################################
###############################################################################
###############################################################################


rule makeJSONForReciprocalBlastPhylogeny:
  input:
    db    = expand("%s/{genome}.blastdb" % (__BLASTDB_OUTDIR__), genome=config["genomes"].keys()),
    prots = [ config["genomes"][genome]["prots"] for genome in config["genomes"].keys() ]
  output:
    json = "%s/recipBlast.json" % __PHYLO_OUTPUT__
  run:
    rbJson = { "args" : {},
               "reciprocalDir" : __RECIPROCALBLAST_OUTPUT__,
               "final" : "%s/reciprocalList.tsv" % __PHYLO_OUTPUT__,
               "evalue" : dconfig["evalue"] }

    referenceBlastDB = "%s/%s.blastdb" % (__BLASTDB_OUTDIR__, config["phylogeny_reference"])
    referenceProts = dconfig["genomes"][dconfig["phylogeny_reference"]]["prots"]
    for genome in dconfig["genomes"]:
      if genome in [ dconfig["phylogeny_reference"], dconfig["querySetName"] ]:
        continue
      #fi
      targetBlastDB = "%s/%s.blastdb" % (__BLASTDB_OUTDIR__, genome)
      targetProts   = dconfig["genomes"][genome]["prots"]

      rbJson["args"]["%s:%s" % (dconfig["phylogeny_reference"], genome)] = { "query" : referenceProts, "db": targetBlastDB }
      rbJson["args"]["%s:%s" % (genome, dconfig["phylogeny_reference"])] = { "query" : targetProts, "db": referenceBlastDB }
    #efor

    with open(output.json, "w") as ofd:
      ofd.write(json.dumps(rbJson, indent=4))
    #ewith

rule runReciprocalBlastPhylogeny:
  input:
    json = rules.makeJSONForReciprocalBlastPhylogeny.output.json
  output:
    reciprocalList = "%s/reciprocalList.tsv" % __PHYLO_OUTPUT__
  threads: 99
  params:
    pc_dir = __PC_DIR__
  shell: """
    {params.pc_dir}/run_reciprocal_blast.sh "{input.json}" {threads} #&> /dev/null
  """

rule makeJSONForClustalo:
  input:
    reciprocalList = rules.runReciprocalBlastPhylogeny.output.reciprocalList
  output:
    json = "%s/clustalo.json"% __PHYLO_OUTPUT__
  run:
    recipBlastFiles = dict([ (r.genomepair, r.file) for r in utils.readColumnFile(input.reciprocalList, "genomepair file")])
    refProts = utils.loadFasta(config["genomes"][dconfig["phylogeny_reference"]]["prots"])

    matches = { genome: {} for genome in dconfig["genomes"] if genome != dconfig["querySetName"] }
    matches[dconfig["phylogeny_reference"]] = { gene : gene for gene in refProts }

    for genome in dconfig["genomes"]:
      if genome in [ dconfig["phylogeny_reference"], dconfig["querySetName"] ]:
        continue
      #fi
      refValt = butils.bestHitPerQuery(butils.readBlastFile(recipBlastFiles["%s:%s" % (config["phylogeny_reference"], genome)], butils.diamondfields))
      altVref = butils.bestHitPerQuery(butils.readBlastFile(recipBlastFiles["%s:%s" % (genome, config["phylogeny_reference"])], butils.diamondfields))

      for seq in refProts:
        bhSeq   = refValt[seq].sseqid if seq in refValt else "__NOMATCH"
        bhrecip = altVref[bhSeq].sseqid if bhSeq in altVref else "__NOMATCH"

        if bhrecip == seq:
          matches[genome][seq] = bhSeq
        #fi
      #efor

      # Find out how many genomes each gene is present in
    geneGenomeCounts = {}
    for gene in refProts:
      geneGenomeCounts[gene] = len([ genome for genome in dconfig["genomes"] if (genome != dconfig["querySetName"]) and (gene in matches[genome]) ])
    #efor

    # Select genes present in at least dconfig["phylogeny_min_prevalence"] % genomes
    selected = []
    for gene in refProts:
      prevalence = (float(geneGenomeCounts[gene]) / float(len(dconfig["genomes"]) -1)) * 100.0
      if prevalence >= (dconfig["phylogeny_min_prevalence"]) :
        selected.append(gene)
      #fi
    #efor

    print("Found %d genes present in atleast %.1f genomes. Will select top %d (based on number of genomes it was found in)." % (len(selected), dconfig["phylogeny_min_prevalence"],  min(len(selected), dconfig["phylogeny_max_genes"]) ))
    print("If you want to change this behaviour, change 'phylogeny_max_genes', and 'phylogeny_min_prevalence' in config.")

    J = { "outdir" : __PHYLO_OUTPUT__,
          "args" : {},
          "final" : rules.alignPhylogeny.output.aln,
          "genomes" : [ g for g in dconfig["genomes"].keys() if g != dconfig["querySetName"] ]
        }

    otherProts = { genome : utils.loadFasta(config["genomes"][genome]["prots"]) for genome in matches }
    for gene in sorted(selected, key=lambda x: geneGenomeCounts[x], reverse=True):
      utils.writeFasta([ (genome, otherProts[genome][matches[genome][gene]]) for genome in sorted(matches.keys()) if (gene in matches[genome]) ], "%s/%s.fasta" % (__PHYLO_OUTPUT__, gene))
      J["args"][gene] = "%s/%s.fasta" % (__PHYLO_OUTPUT__, gene)
    #efor

    with open(output.json, "w") as ofd:
      ofd.write(json.dumps(J, indent=4))
    #ewith
          
rule alignPhylogeny:
  input:
    json = rules.makeJSONForClustalo.output.json
  output:
    aln = "%s/allAligned.fasta" % __PHYLO_OUTPUT__
  threads: 99
  params:
    pc_dir = __PC_DIR__
  shell: """
    {params.pc_dir}/run_clustalo_aligner.sh "{input.json}" {threads} #&> /dev/null
  """

rule fastTreePhylogeny:
  input:
    aln = rules.alignPhylogeny.output.aln
  output:
    newick = "%s/phylogeny.newick" % __PHYLO_OUTPUT__
  threads: 99
  params:
    fasttreeParams = dconfig["fasttree_params"]
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    export OMP_NUM_THREADS={threads}
    FastTree {params.fasttreeParams} < "{input.aln}" > {output.newick}
  """

rule rerootedPhylogeny:
  input:
    newick = rules.fastTreePhylogeny.output.newick
  output:
    newick = "%s/pyhlogeny.rerooted.newick" % __PHYLO_OUTPUT__
  params:
    outgroups = dconfig["phylogeny_outgroups"]
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    nw_reroot "{input.newick}" {params.outgroups} > "{output.newick}"
  """
rule phylogeny: 
  input:
    newick = rules.rerootedPhylogeny.output.newick

###############################################################################
###############################################################################
###############################################################################

#     1  149_Ascsa1:149_Ascsa1|4785
#     2  e03ea43b0cfb2e7ff623bc050d3ac750
#     3  194
#     4  Pfam
#     5  PF00071
#     6  Ras family
#     7  7
#     8  177
#     9  8.4E-53
#    10  T
#    11  15-09-2017
#    12  IPR001806
#    13  Small GTPase superfamily

rule interproscan_initial:
  input:
    final = lambda wildcards: dconfig["queries"][wildcards.query]
  output:
    annot = "%s/ipr.initial.{query}.tsv" % __VALIDATION_OUTDIR__
  shell: """
    cat {input.final} \
     | tr -d '*' \
     > {output.annot}.input
    interproscan.sh -appl Pfam -f TSV --iprlookup -i {output.annot}.input -o {output.annot}
  """

rule interproscan_final:
  input:
    final = lambda wildcards: "%s/%s.fasta" % (__FINAL_OUTDIR__, wildcards.query)
  output:
    annot = "%s/ipr.final.{query}.tsv" % __VALIDATION_OUTDIR__
  shell: """
    cat {input.final} \
     | tr -d '*' \
     > {output.annot}.input
    interproscan.sh -appl Pfam -f TSV --iprlookup -i {output.annot}.input -o {output.annot}
  """

rule validationScorePerGene:
  input:
    iAnnots = expand("%s/ipr.initial.{query}.tsv" % __VALIDATION_OUTDIR__, query=sorted(dconfig["queries"].keys())),
    fAnnots = expand("%s/ipr.final.{query}.tsv" % __VALIDATION_OUTDIR__, query=sorted(dconfig["queries"].keys())),
    final  = expand("%s/{query}.fasta" % __FINAL_OUTDIR__, query=sorted(dconfig["queries"].keys()))
  output:
    scores = "%s/geneScores.tsv" % __VALIDATION_OUTDIR__
  run:
    C = {}
    for query, iAnnotFile, fAnnotFile in zip(sorted(dconfig["queries"].keys()), input.iAnnots, input.fAnnots):
      initialAnnots = set([ a.signatureid for a in utils.readColumnFile(iAnnotFile, "seqid hash seqlen analysis signatureid signaturedesc start stop score status date iprid iprdesc") ])
      finalAnnots = utils.readColumnFile(fAnnotFile, "seqid hash seqlen analysis signatureid signaturedesc start stop score status date iprid iprdesc")

      for a in finalAnnots:
        if (query, a.seqid) not in C:
          C[(query, a.seqid)] = [ initialAnnots, set([]) ]
        #fi
        initial, final = C[(query, a.seqid)]
        final.add(a.signatureid)

      #efor
    #efor

    for query, seqFile in zip(sorted(dconfig["queries"].keys()), input.final):
      seqs = utils.loadFasta(seqFile)
      for seq in seqs:
        if (query, seq) not in C:
          C[(query, seq)] = [ set([]), set([]) ]
        #fi
      #efor
    #efor
        
 
    with open(output.scores, "w") as ofd:
      for (query, geneid) in C:
        initial, final = C[(query, geneid)]
        nInitial = len(initial)
        nOverlap = len(initial & final)
        ofd.write("%s\t%s\t%d\t%d\t%d\t%f\n" % (query, geneid, nInitial, len(final), nOverlap, ( nOverlap / float(nInitial) ) if nInitial > 0 else 1) )
      #efor
    #ewith

rule validationScores:
  input:
    iAnnots = expand("%s/ipr.initial.{query}.tsv" % __VALIDATION_OUTDIR__, query=sorted(dconfig["queries"].keys())),
    fAnnots = expand("%s/ipr.final.{query}.tsv" % __VALIDATION_OUTDIR__, query=sorted(dconfig["queries"].keys())),
    final  = expand("%s/{query}.fasta" % __FINAL_OUTDIR__, query=sorted(dconfig["queries"].keys()))
  output:
    valTable = "%s/validation.tsv" % __VALIDATION_OUTDIR__,
    valTableTrans = "%s/validation.trans.tsv" % __VALIDATION_OUTDIR__
  run:
    import statistics
    S = { "NGenes" : "Number of genes", "NAnnots" : "Number of annotations", "Average" : "Average", "Maximum" : "Maximum", "Minimum" : "Minimum", "Median" : "Median" }
    A = { query : { seqid : [] for seqid in utils.loadFasta(final) } for (query, final) in zip(sorted(dconfig["queries"]), input.final) }
    for query, iAnnotFile, fAnnotFile in zip(sorted(dconfig["queries"].keys()), input.iAnnots, input.fAnnots):
      initialAnnot = utils.readColumnFile(iAnnotFile, "seqid hash seqlen analysis signatureid signaturedesc start stop score status date iprid iprdesc")
      finalAnnot = utils.readColumnFile(fAnnotFile, "seqid hash seqlen analysis signatureid signaturedesc start stop score status date iprid iprdesc")
      relAnnots = set([ a.signatureid for a in initialAnnot ])
      for a in finalAnnot:
        if a.signatureid not in relAnnots:
          continue
        #fi
        A[query][a.seqid].append(a.signatureid)
        if a.signatureid not in S:
          S[a.signatureid] = a.signaturedesc
        #fi
      #efor
    #efor

    V = {}
    for query in A:
      nseq = float(len(A[query].keys()))
      counts = []
      for sig in S:
        nSigQA = sum([ 1 if sig in A[query][seqid] else 0 for seqid in A[query] ])
        V[(query, sig)] = -1 if nSigQA == 0 else 100*(nSigQA / nseq)
        if nSigQA > 0:
          counts.append(V[(query, sig)])
        #fi
      #efor
      V[(query, "NGenes")] = len(A[query].keys())
      V[(query, "NAnnots")] = len(counts)
      V[(query, "Average")] = statistics.mean(counts) if len(counts) > 0 else 100
      V[(query, "Maximum")] = max(counts) if len(counts) > 0 else 100
      V[(query, "Minimum")] = min(counts) if len(counts) > 0 else 100
      V[(query, "Median")] = statistics.median(counts) if len(counts) > 0 else 100
    #efor

    with open(output.valTable, "w") as ofd:
      ofd.write("Annotation\tDescription\t%s\n" % "\t".join(sorted(dconfig["queries"].keys())))
      for annot in S.keys():
        ofd.write("%s\t%s\t%s\n" % (annot, S[annot], "\t".join([ "%.0f" % V[(query, annot)] for query in sorted(dconfig["queries"].keys()) ]) ))
      #efor
    #ewith

    with open(output.valTableTrans, "w") as ofd:
      keyS = S.keys()
      ofd.write("Query\t%s\n" % "\t".join(keyS))
      ofd.write("\t%s\n" % "\t".join([ S[s] for s in keyS]))
      for query in sorted(dconfig["queries"].keys()):
        ofd.write("%s\t%s\n" % (query, "\t".join([ "%.0f" % V[(query, annot)] for annot in keyS])) )
      #efor
    #ewith

rule filterBasedOnValidation:
  input:
    scores = rules.validationScorePerGene.output.scores,
    final = rules.final.input.final,
    cmpTable = rules.compareToInitial.output.cmpTable,
    nodesCSV = rules.matchedListToGraph.output.nodesCSV
  output:
    cmpTable = "%s/validated.cmpTable.tsv"% __VALIDATION_OUTDIR__,
    nodesCSV = "%s/validated.nodes.tsv" % __VALIDATION_OUTDIR__,
    removed  = "%s/validated.removed.tsv" % __VALIDATION_OUTDIR__
  run:
    ml = utils.indexListBy(utils.readColumnFile("%s/allMatched.tsv" % __ITERATION_OUTDIR__, "query iteration origin target keep", types="str int str str int"), lambda x: x.iteration)
    vScore = utils.indexListBy(utils.readColumnFile(input.scores, "query geneid ninitial nfinal noverlap score", types="str str int int int float"), lambda x: x.geneid)

    selected = { m.origin: m.query for m in ml[0] }
    removed = {}
    for iteration in range(max([ int(k) for k in ml.keys()]) + 1):
      for m in ml[iteration]:
        if m.keep < 1:
          continue
        elif m.target not in vScore:
          selected[m.target] = m.query
        elif (dconfig["validation_filter_strategy"] == "min_score") and (vScore[m.target][0].score > dconfig["validation_filter_min_score"]):
          selected[m.target] = m.query
        elif m.origin not in selected:
          removed[m.query] = removed.get(m.query,[]) + [ m.target ]
          continue
        elif (dconfig["validation_filter_strategy"] == "min_score_along_graph") and (vScore[m.target][0].score > dconfig["validation_filter_min_score"]):
          selected[m.target] = m.query
        else:
          removed[m.query] = removed.get(m.query,[]) + [ m.target ]
        #fi
      #efor
    #efor

    with open(output.removed, "w") as ofd:
      for query in removed:
        for gene in removed[query]:
          ofd.write("%s\t%s\n"% (query, gene))
        #efor
      #efor
    #ewith

    querySpecies = { query: set([]) for query in dconfig["queries"] }
    for s in selected:
      genome, gene = utils.splitgg(s)
      querySpecies[selected[s]].add(genome)
    #efor
      

    queries = sorted(dconfig["queries"].keys())
    with open(input.cmpTable, "r") as ifd:
      with open(output.cmpTable, "w") as ofd:
        for line in ifd:
          line = line.rstrip().split('\t')
          if line[0] == 'Genome':
            ofd.write('\t'.join(line) + '\n')
          else:
            genome = line[0]
            newLine = [ line[0] ]
            for query, qual in zip(queries, line[1:]):
              if genome not in querySpecies[query]:
                if qual == "-1":
                  newLine.append(qual)
                else:
                  newLine.append('-' + qual)
                #fi
                #print("I remove %s from query %s" % (genome, query))
              else:
                newLine.append(qual)
              #fi
            #efor
            ofd.write('\t'.join(newLine) + '\n')
          #fi
        #efor
      #ewith
    #ewith

    nodes = utils.readColumnFile(input.nodesCSV, "id query label iteration", delimiter = ',', types="str str str int", skip=1)
    with open(output.nodesCSV, "w") as ofd:
      ofd.write("Id,query,label,iteration,removed,vscore\n")
      for node in nodes:
        rm = 0
        if node.label in removed.get(node.query, []):
          rm = 1
        #fi
        ofd.write("%s,%s,%s,%d,%d,%f\n" % (node.id, node.query, node.label, node.iteration, rm, vScore[node.label][0].score if node.label in vScore else 0.0))
      #efor
    #ewith

rule itolAnnotPerQuery:
  input:
    cmpTable = rules.filterBasedOnValidation.output.cmpTable,
    template = "%s/raw_itol_colored_gradients.txt" % __PC_DIR__
  output:
    itol = "%s/validated_itol/itol.{query}.txt" % __PHYLO_OUTPUT__
  run:
    queries = sorted(dconfig["queries"].keys())
    idx = queries.index(wildcards.query)
    data = []
    with open(input.cmpTable, "r") as ifd:
      for line in ifd:
        line = line.rstrip().split('\t')
        if line[0] == 'Genome':
          continue
        else:
          genome = line[0]
          qscore = int(line[1+idx])
          data.append("%s %d" % (genome, 0 if qscore < 0 else qscore))
        #fi
      #efor
    #ewith

    with open(input.template, "r") as tfd:
      tContent = tfd.read()
      with open(output.itol, "w") as ofd:
        ofd.write(tContent % (wildcards.query, "\n".join(data)))
      #ewith
    #ewith
    #

rule itolAnnot:
  input:
    annot = expand("%s/validated_itol/itol.{query}.txt"% __PHYLO_OUTPUT__, query=dconfig["queries"].keys())

rule validate:
  input:
    scores = rules.validationScores.output.valTable,
    filt   = rules.filterBasedOnValidation.output.cmpTable,
    itolannot = rules.itolAnnot.input.annot
