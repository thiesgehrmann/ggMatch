
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

###############################################################################

querySetProts = "%s/querySet.fasta" % __ITERATION_OUTDIR__

dconfig["genomes"][dconfig["querySetName"]] = {"prots" : querySetProts }

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
      {params.pc_dir}/run_iteration.sh {input.json} {threads}
    fi
  """
    
#rule runAllIterations:
#  input:
#    groups = lambda wildcards: expand("%s/%s/iteration_{iteration}.fasta" % (__ITERATION_OUTDIR__, wildcards.query), iteration=[ str(x) for x in range(1,10)])

rule moveFinal:
  input:
    groups = lambda wildcards: expand("%s/%s/iteration_{iteration}.fasta" % (__ITERATION_OUTDIR__, wildcards.query), iteration=[ str(x) for x in range(1,10)]),
    lastIter = lambda wildcards: "%s/%s/iteration_9.fasta" % (__ITERATION_OUTDIR__, wildcards.query),
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

rule interproScan:
  input:
    final = lambda wildcards: "%s/%s.fasta" % (__FINAL_OUTDIR__, wildcards.query)
  output:
    annot = "%s/ipr.{query}.tsv" % __VALIDATION_OUTDIR__
  shell: """
    cat {input.final} \
     | tr -d '*' \
     > {output.annot}.input
    interproscan.sh -appl Pfam -f TSV --goterms --iprlookup -i {output.annot}.input -o {output.annot}
  """

rule validationScores:

rule validate:
  input:
    scores = rules.validationScores.output.scores
