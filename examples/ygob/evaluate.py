#!/usr/bin/env python
import matplotlib
matplotlib.use('agg')

import ygob as ygob

import numpy as np
import sys
import json
import os.path
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

J = json.load(open(sys.argv[1], "r"))
P = ygob.readPillars(sys.argv[2])

results = []

def evaluate(genes, pillar, removed=set([])):

  overlap = [ gene for (genome, gene) in [ ygob.splitgg(g) for g in genes ] if gene in pillar[genome] ]

  precision = len(overlap) / float(len(genes))
  recall    = len(overlap) / float(len([ k for k in pillar if len(pillar[k]) > 0]))

  f1 = 2 * (precision * recall) / (precision + recall)

  #print(len(overlap), precision, recall, f1)

  return (precision, recall, f1)

removed = []
if os.path.isfile("output/run/validation/validated.removed.tsv"):
  with open("output/run/validation/validated.removed.tsv", "r") as ifd:
    for line in ifd:
      removed.append(line.split("\t")[1])
    #efor
  #ewith
#fi

D = []

for query in J["queries"]:
  finalFile = "output/run/final/%s.fasta" % query
  queryID = int(query.split('_')[1])
  queryPillar = P[queryID]

  F = ygob.loadFasta(finalFile)


  evaluation = evaluate(F.keys(), queryPillar)
  D.append((query, "Precision", evaluation[0]))
  D.append((query, "Recall",    evaluation[1]))
  D.append((query, "F1",        evaluation[2]))
  results.append(evaluation)
  print(query, evaluation)
#efor

D = np.array(D)

results = np.matrix(results)

print("Average: ", np.mean(results, axis=0))
print("St dev:  ", np.std(results, axis=0))


DF = pd.DataFrame(D, columns = ["query", "measure", "value" ])
DF['value'] = DF['value'].apply(pd.to_numeric)
print(DF)
ax = sns.violinplot(x="measure", y="value", data=DF, cut=0)
plt.savefig("evaluation.png")
plt.savefig("evaluation.svg")
