
# Config has format
# {
#   blastdbs: { genome1 : "blastDBFile" },
#   prots: { genome1 : "genome1prots.fasta" },
#   queries: "queryFile",
#   outdir: "outdir",
#   final: "finalFile"
# }

import utils as utils
import blastutils as butils

###############################################################################

__OUTDIR__ = config["outdir"]

###############################################################################

rule all:
  input:
    file = config["final"]

###############################################################################

rule runBlast:
  input:
    db = lambda wildcards: config["blastdbs"][wildcards.genome],
    query = config["queries"]
  output:
    res = "%s/blastres.{genome}.tsv" % __OUTDIR__
  params:
    blastfields = butils.blastfields
  conda: "conda.yaml"
  shell: """
    blastp -query {input.query} -db {input.db} -outfmt "6 {params.blastfields}" -out {output.res} -num_threads 1 -max_target_seqs 1 -evalue 1e-50
  """

###############################################################################

rule iterationFindBest:
  input:
    res   = lambda wildcards: expand("%s/blastres.{genome}.tsv" % (__OUTDIR__), genome=sorted(config["prots"].keys())),
    prots = [ config["prots"][genome] for genome in sorted(config["prots"].keys()) ]
  output:
    fas = "%s/newOrthologs.fasta" % (__OUTDIR__)
  run:
    F = utils.loadFasta(config["queries"])
    for (genome, resFile, protFile) in zip(sorted(config["blastdbs"].keys()), input.res, input.prots):
      bRes = butils.readBlastFile(resFile, butils.blastfields)

      if len(bRes) == 0:
        continue
      #fi

      prots = utils.loadFasta(protFile)

      #Initially just take the best hit, if there is any
      best = butils.bestHit(bRes)
      F["%s:%s" % (genome, best.sseqid)] = prots[best.sseqid]
    #efor
    print(F)
    utils.writeFasta(F.items(), output.fas)

###############################################################################

rule final:
  input:
    file = rules.iterationFindBest.output.fas
  output:
    file = config["final"]
  shell: """
    mv {input.file} {output.file}
  """

###############################################################################
