import blastutils as butils

rule all:
  input:
    final = config["final"]

###############################################################################

rule reciprocalBlast:
  input:
    query = lambda wildcards: config["args"][wildcards.pair]["query"],
    db    = lambda wildcards: config["args"][wildcards.pair]["db"]
  output:
    res = "%s/blast.{pair}.tsv" % config["reciprocalDir"]
  threads: 4
  conda: "conda.yaml"
  params:
    blastfields = butils.blastfields,
    evalue = config["evalue"]
  shell: """
    #blastp -query {input.query} -db {input.db} -outfmt "6 {params.blastfields}" -out {output.res} -num_threads {threads} -max_target_seqs 1
    diamond blastp -d {input.db} -q {input.query} -o {output.res} --max-target-seqs 1 --evalue {params.evalue}
  """

###############################################################################

rule final:
  input:
    res = expand("%s/blast.{pair}.tsv" % config["reciprocalDir"], pair=sorted(config["args"].keys()))
  output:
    final = config["final"]
  run:
    with open(output.final, "w") as ofd:
      for pair in config["args"].keys():
        res = "%s/blast.%s.tsv" % (config["reciprocalDir"], pair)
        ofd.write("%s\t%s\n" % (pair, res))
      #efor
    #ewith
