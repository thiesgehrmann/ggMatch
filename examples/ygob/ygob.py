#!/usr/bin/env python

from collections import namedtuple
import gzip
import json
import sys

###############################################################################

def loadFasta(fastaFile):

  F = {'': 0}

  current_seq = ""
  with (gzip.open(fastaFile, "r") if fastaFile[-2:] == "gz" else open(fastaFile, "r")) as fd:
    for line in fd:
      line = line.strip()
      if len(line) == 0:
        continue
      #fi
      if line[0] == '>':
        current_seq = line[1:].split(' ')[0]
        F[current_seq] = ""
      else:
        F[current_seq] = F[current_seq] + line.strip()
      #fi
  #ewith
  F.pop("", None)
  return F
#edef

###############################################################################

def writeFasta(fasta, outFile, linelength=80):
  with open(outFile, "w") as ofd:
    for  (name, sequence) in fasta:
      ofd.write(">%s\n" % name)
      ofd.write("%s\n" % '\n'.join([sequence[i:i+linelength] for i in range(0, len(sequence), linelength)]))
    #efor
  #ewith
#edef

###############################################################################

def readColumnFile(filename, columnNames, delimiter='\t', types=""):
  import csv
  L = []
  typeFunctions = { "str" : lambda x: str(x),
                    "int" : lambda x: int(x),
                    "float" : lambda x: float(x) }

  if types != "":
    types = [ typeFunctions[c] for c in types.split(" ") ]
  #fi

  lineType = namedtuple("lineType", columnNames)
  with open(filename, "r") as ifd:
    reader = csv.reader(ifd, delimiter=delimiter)
    for row in reader:
      if row[0][0] == '#':
        continue
      #fi
      if len(types) == len(row):
        row = [ tf(v) for (tf, v) in zip(types, row) ]
      #fi
      L.append(lineType(*row[:len(columnNames)]))
    #efor
  #ewith
  return L
#edef

###############################################################################

def splitgg(gg):
  ggs = gg.split(':')

  genome = ggs[0]
  gene   = ":".join(ggs[1:])

  if len(ggs) == 1:
    return (None, gene)
  #fi

  return (genome, gene)

###############################################################################

genomes = [
  ("lwalt", "Kwal_", "L. waltii"),
  ("egoss", "A",     "E. gossypii"),
  ("cglab", "CAGL0", "C. glabrata"),
  ("ecymb", "Ecym_", "E. cymbalariae"),
  ("klact", "KLLA0", "K. lactis"),
  ("lther", "KLTH0", "L. thermotolerans"),
  ("knaga", "KNAG0", "K. naganishii"),
  ("vpoly", "Kpol_", "V. polyspora"),
  ("ncast", "NCAS0", "N. castellii"),
  ("ndair", "NDAI0", "N. dairenensis"),
  ("lkluy", "SAKL0", "L. kluyveri"),
  ("skudr", "Skud_", "S. kudriavzevii"),
  ("smika", "Smik_", "S. mikatae"),
  ("sbaya", "Suva_", "S. bayanus var. uvarum"),
  ("tblat", "TBLA0", "T. blattae"),
  ("tdelb", "TDEL0", "T. delbrueckii"),
  ("tphaf", "TPHA0", "T. phaffii"),
  ("zroux", "ZYRO0", "Z. rouxii"),
  ("kafri", "KAFR0", "K. africana"),
  ("scere", "Y",     "S. cerevisiae") ]

pillarGenomes = [
  ("vpoly", "V. polyspora Position 1"),
  ("tphaf", "T. phaffii Position 1"),
  ("tblat", "T. blattae Position 1"),
  ("ndair", "N. dairenensis Position 1"),
  ("ncast", "N. castellii Position 1"),
  ("knaga", "K. naganishii Position 1"),
  ("kafri", "K. africana Position 1"),
  ("cglab", "C. glabrata Position 1"),
  ("sbaya", "S. bayanus var. uvarum Position 1"),
  ("skudr", "S. kudriavzevii Position 1"),
  ("smika", "S. mikatae Position 1"),
  ("scere", "S. cerevisiae Position 1"),
  ("ago",   "Ancestral Gene Order Position 1"),
  ("zroux", "Z. rouxii Position 1"),
  ("tdelb", "T. delbrueckii Position 1"),
  ("klact", "K. lactis Position 1"),
  ("egoss", "E. gossypii  Position 1"),
  ("ecymb", "E. cymbalariae Position 1"),
  ("lkluy", "L. kluyveri Position 1"),
  ("lther", "L. thermotolerans Position 1"),
  ("lwalt", "L. waltii Position 1"),
  ("scere", "S. cerevisiae Position 2"),
  ("smika", "S. mikatae Position 2"),
  ("skudr", "S. kudriavzevii Position 2"),
  ("sbaya", "S. bayanus var. uvarum Position 2"),
  ("cglab", "C. glabrata Position 2"),
  ("kafri", "K. africana Position 2"),
  ("knaga", "K. naganishii Position 2"),
  ("ncast", "N. castellii Position 2"),
  ("ndair", "N. dairenensis Position 2"),
  ("tblat", "T. blattae Position 2"),
  ("tphaf", "T. phaffii Position 2"),
  ("vpoly", "V. polyspora Position 2")
]

###############################################################################

def readPillars(pillarsFile):
  P = []
  with open(pillarsFile, "r") as ifd:
    for line in ifd:
      p = { g[0]: [] for g in genomes }
      genes = line.split('\t')
      for (genome, name), gene in zip(pillarGenomes, genes):
        if (genome in p) and ('--' not in gene):
          p[genome].append(gene)
        #fi
      #efor
      P.append(p)
    #efor
  #ewith
  return P
#edef

###############################################################################

if __name__ == "__main__":

  fastaFile = sys.argv[1]
  pillarsFile = sys.argv[2]
  
  F = loadFasta(fastaFile)
  J = { "genomes" : {},
        "queries" : {},
        "outdir" : "./examples/ygob/output" }

  # Split the protein sequence file into different files
  for genome, prefix, name in genomes:
    keys = [ sid for sid in F.keys() if sid.startswith(prefix) ]
    writeFasta([ (sid.split(' ')[0], F[sid]) for sid in keys ], "genomes/%s.fasta" % genome)
    J["genomes"][genome] = { "prots" : "examples/ygob/genomes/%s.fasta" % genome,
                             "name" : name }
  #efor
  
  
  P = readPillars(pillarsFile)
  
  totali = 0
  for pi, p in enumerate(P):
    if len(p["scere"]) > 0:
      gene = p["scere"][0]
      writeFasta([ ("scere:%s" % gene, F[gene])], "queries/query_%d.fasta" % pi)
      J["queries"]["query_%d" % pi] = "examples/ygob/queries/query_%d.fasta" % pi
      totali += 1
    #fi
  
    if totali >= 40:
      break
  #done
  
  
  
  json.dump(J, open("config.json", "w"), indent=4)
