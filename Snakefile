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
__DOT_OUTDIR__ = "%s/dot" % __RUN_DIR__
###############################################################################

rule makeBlastDB:
  input:
    proteins = lambda wildcards: dconfig["genomes"][wildcards.genome]["prots"]
  output:
    db = "%s/{genome}.blastdb" % __BLASTDB_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    makeblastdb -dbtype prot -in {input.proteins} -out {output.db}
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
          "queryName" : wildcards.query,
          "iteration" : int(wildcards.iteration),
          "queries" : input.queries,
          "outdir" : "%s/%s/iteration_%d.run" % ( __ITERATION_OUTDIR__, query, iteration),
          "final" : "%s/%s/iteration_%d.fasta" % ( __ITERATION_OUTDIR__, query, iteration),
          "matchedList" : "%s/allMatched.tsv" % (__ITERATION_OUTDIR__),
          "needNext" : needNext }
    with open(output.json, "w") as ofd:
      ofd.write(json.dumps(C, indent=4))

rule initIteration0:
  input:
    query = lambda wildcards: dconfig["queries"][wildcards.query]
  output:
    query = "%s/{query}/iteration_0.fasta" % (__ITERATION_OUTDIR__)
  shell: """
     cat "{input.query}" \
      | sed -e 's/^>/>__INPUTQUERY__:/' \
      > "{output.query}"
  """

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
      {params.pc_dir}/run_iteration.sh {input.json} {threads} &> /dev/null
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

rule matchedListToDot:
  input:
    files = rules.final.input.final
  output:
    dot = "%s/matchedList.dot" % __DOT_OUTDIR__
  run:
    ml = utils.readColumnFile("%s/allMatched.tsv" % __ITERATION_OUTDIR__, "query iteration queryMatch genome gene", types="str int str str str")

    visited = set([])
    colors = dict(zip(range(0,10), ["gray", "black", "red", "salmon", "brown", "orange", "yellow", "darkgreen", "limegreen", "skyblue" ]))

    with open(output.dot, "w") as ofd:
      ofd.write("digraph G{\n")
      ofd.write("  subgraph clusteregend{\n")
      #ofd.write("    node [style=filled, color=white]\n")
      ofd.write("    style = filled;\n")
      ofd.write("    color = lightgrey;\n")
      ofd.write("    label = \"Legend\";\n")
      for iteration in range(0,10):
        if iteration == 0:
          ofd.write("    iteration0 [color=red,peripheries=1];\n")
        else:
          ofd.write("    iteration%d [peripheries=%d];\n" % (iteration, iteration))
        #fi
      #efor

      for iteration in range(1,10):
        ofd.write("    iteration%d -> iteration%d;\n" % (iteration-1, iteration))
      #efor
      ofd.write("  }\n\n")

      nodes = {}
      nodeAttr = {}
      edges = {}
      nodeIDcounter = 1
      for match in ml:
        originGenome = match.queryMatch.split(":")[0]
        originGene   = ':'.join(match.queryMatch.split(":")[1:])
        targetGenome = match.genome
        targetGene   = match.gene

        query = match.query
        iteration = match.iteration
        origin = "%s:%s" % (originGenome, originGene)
        target = "%s:%s" % (targetGenome, targetGene)

        print(origin, target)

        if query not in nodes:
          nodes[query] = {}
          nodeAttr[query] = {}
          edges[query] = []
        #fi
        if origin not in nodes[query]:
          nodes[query][origin] = "node%d" % nodeIDcounter
          nodeIDcounter += 1
          if iteration == 1:
            nodeAttr[query][origin] = "[peripheries=1,color=red]"
          else:
            nodeAttr[query][origin] = "[peripheries=%d]" % iteration
          #fi
        #fi

        if target not in nodes[query]:
          nodes[query][target] = "node%d" % nodeIDcounter
          nodeIDcounter += 1
          nodeAttr[query][target] = "[peripheries=%d]" % iteration
        #fi

        edges[query].append( (nodes[query][origin], nodes[query][target]))

      #efor

      for query in nodes:
        ofd.write("  subgraph cluster%s {\n" % query)
        ofd.write("    style = filled;\n")
        ofd.write("    color = lightgrey;\n")
        ofd.write("    label = \"%s\";\n" % query)
        for node in nodes[query]:
          ofd.write("    %s %s;\n" % ( nodes[query][node], nodeAttr[query][node]))
        #efor
        for (origin, target) in edges[query]:
          ofd.write("    %s -> %s;\n" % (origin, target) )
        #efor
        ofd.write("  }\n")
      #efor

      ofd.write("}")

rule matchedDotToPdf:
  input:
    dot = rules.matchedListToDot.output.dot,
  output:
    pdf = "%s/matchedList.pdf" % __DOT_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    dot -Tpdf "{input.dot}" > "{output.pdf}"
  """
###############################################################################
