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

###############################################################################

noGenomeProts = "%s/noGenome.prots.fasta" % __ITERATION_OUTDIR__

dconfig["genomes"]["NoGenome"] = {"prots" : noGenomeProts }

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
          "reciprocalDir" : "%s/reciprocalBlasts" % __ITERATION_OUTDIR__,
          "outdir" : "%s/%s/iteration_%d.run" % ( __ITERATION_OUTDIR__, query, iteration),
          "final" : "%s/%s/iteration_%d.fasta" % ( __ITERATION_OUTDIR__, query, iteration),
          "matchedList" : "%s/allMatched.tsv" % (__ITERATION_OUTDIR__),
          "evalue" : dconfig["evalue"],
          "needNext" : needNext,
          "pcDir" : __PC_DIR__ }
    with open(output.json, "w") as ofd:
      ofd.write(json.dumps(C, indent=4))

rule genNoGenomeProts:
  input:
    files = [ dconfig["queries"][query] for query in dconfig["queries"].keys() ]
  output:
    noGenomeProts = dconfig["genomes"]["NoGenome"]["prots"]
  shell: """
    cat {input.files} > {output.noGenomeProts}
  """

rule initIteration0:
  input:
    query = lambda wildcards: dconfig["queries"][wildcards.query],
    noGenomeProts = rules.genNoGenomeProts.output.noGenomeProts
  output:
    query = "%s/{query}/iteration_0.fasta" % (__ITERATION_OUTDIR__)
  run:
    F  = utils.loadFasta(input.query)
    NF = {}
    for ident in F:
      genome = ident.split(':')[0]
      if genome not in dconfig["genomes"]:
        NF["NoGenome:%s" % ident] = F[ident]
      else:
        NF[ident] = F[ident]
      #fi
    #efor
    utils.writeFasta(NF.items(), output.query)

    with open("%s/allMatched.tsv" % __ITERATION_OUTDIR__, "a+") as ofd:
      for seed in NF:
        ofd.write("%s\t%d\tQUERY:%s\t%s\t%d\n" % (wildcards.query, 0, wildcards.query, seed, 1))
      

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
    genomes = [ g for g in dconfig["genomes"] if g != "NoGenome" ]
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
          
          

      
    
