#!/usr/bin/env python

import ygob as ygob

import numpy as np
import sys
import json

J = json.load(open(sys.argv[1], "r"))
P = ygob.readPillars(sys.argv[2])

results = []

def evaluate(genes, pillar):

  overlap = [ gene for (genome, gene) in [ ygob.splitgg(g) for g in genes ] if gene in pillar[genome] ]

  print(genes)
  print(overlap)
  print(pillar)

  precision = len(overlap) / float(len(genes))
  recall    = len(overlap) / float(len([ k for k in pillar if len(pillar[k]) > 0]))

  f1 = 2 * (precision * recall) / (precision + recall)

  #print(len(overlap), precision, recall, f1)

  return (precision, recall, f1)


for query in J["queries"]:
  finalFile = "output/run/final/%s.fasta" % query
  queryID = int(query.split('_')[1])
  queryPillar = P[queryID]

  F = ygob.loadFasta(finalFile)


  evaluation = evaluate(F.keys(), queryPillar)
  results.append(evaluation)
  print(query, evaluation)
#efor

results = np.matrix(results)

print("Average: ", np.mean(results, axis=0))
print("St dev:  ", np.std(results, axis=0))

