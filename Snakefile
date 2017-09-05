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
          "queries" : input.queries,
          "outdir" : "%s/%s/iteration_%d.run" % ( __ITERATION_OUTDIR__, query, iteration),
          "final" : "%s/%s/iteration_%d.fasta" % ( __ITERATION_OUTDIR__, query, iteration),
          "needNext" : needNext }
    with open(output.json, "w") as ofd:
      ofd.write(json.dumps(C, indent=4))

rule initIteration0:
  input:
    query = lambda wildcards: dconfig["queries"][wildcards.query]
  output:
    query = "%s/{query}/iteration_0.fasta" % (__ITERATION_OUTDIR__)
  shell: """
     ls "{input.query}"
     cp "{input.query}" "{output.query}"
  """

rule runIteration:
  input:
    blastdbs = rules.blastDBs.input.dbs,
    prevQuery = lambda wildcards: "%s/%s/iteration_%s.fasta" % (__ITERATION_OUTDIR__, wildcards.query, int(wildcards.iteration)-1),
    json = lambda wildcards: "%s/%s/%s.json" % (__ITERATION_OUTDIR__, wildcards.query, int(wildcards.iteration))
  output:
    nextQuery = "%s/{query}/iteration_{iteration,[1-9]}.fasta" % __ITERATION_OUTDIR__
  threads: 10
  params:
    pc_dir = __PC_DIR__
  shell: """
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
    lastIter = lambda wildcards: "%s/%s/iteration_9.fasta" % (__ITERATION_OUTDIR__, wildcards.query)
  output:
    finalLoc = "%s/{query}.fasta" % (__FINAL_OUTDIR__)
  shell: """
    ln -s "{input.lastIter}" "{output.finalLoc}"
  """

rule final:
  input:
    final = expand("%s/{query}.fasta" % __FINAL_OUTDIR__, query=dconfig["queries"].keys())

