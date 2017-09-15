import utils as utils
import blastutils as butils


import json

###############################################################################

__OUTDIR__ = config["outdir"]
__ALIGNMENT_OUTDIR__ = "%s/alignment" % __OUTDIR__

###############################################################################

rule all:
  input:
    file = config["final"]

rule clustalo:
  input:
    fasta = lambda wildcards: config["args"][wildcards.gene]
  output:
    aln = "%s/aln.{gene}.fasta" % __ALIGNMENT_OUTDIR__
  conda: "conda.yaml"
  threads: lambda x: int(max(2, 99 / len(config["args"].keys())))
  shell: """
    clustalo -i "{input.fasta}" --threads 1 -o "{output.aln}"
  """

rule mergeAlignments:
  input:
    alignments = expand("%s/aln.{gene}.fasta" % __ALIGNMENT_OUTDIR__, gene=config["args"].keys())
  output:
    final = config["final"]
  run:
    final = { genome : "" for genome in config["genomes"] }
    for alnFile in input.alignments:
      aln = utils.loadFasta(alnFile)
      alnLen = len(aln[list(aln.keys())[0]])
      for genome in sorted(config["genomes"]):
        final[genome] += aln.get(genome, '-' * alnLen)
      #efor
    #efor

    utils.writeFasta(final.items(), output.final)
  
