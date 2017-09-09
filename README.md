# ggMatch
Greedy Gene Matching tool.

ggMatch finds reciprocal best blast/diamond hits across a large number of genomes, in an iterative fashion. 
The iterative nature removes the need to do the dreaded all-vs-all blast between all genomes.
The figure below shows an example gene matching process for a single gene against 298 fungal genomes extracted from the JGI.
The black node represents the original query sequence, and edges between other nodes represent high quality reciprocal best blast hits.
Each node color represents a different iteration.
In the first iteration, we discover the yellow genes, and in later iterations we discover other genes based on the newly discovered genes.
A traditional reciprocal best hit search would have revealed only the yellow nodes in this network.

![Example gene graph created by ggMatch](images/process.png)

## Description of method

COMING SOON

## Usage

### Dependencies

 - Snakemake
 - Conda

### Getting started

You can download and run the example dataset with the following commands:

```bash
  git clone https://github.com/thiesgehrmann/ggMatch.git
  cd ggMatch
  ./ggMatch example/config.json
```

### Output

The tool outputs three files:
 * *outdir/run/compare/cmpTable.tsv* : A matrix of similarity scores to the original query for each species (-1 if absent)

  This score is determined by the number of positive scoring matches in the alignment divided by the length of the query

 * *outdir/run/graph/nodes.tsv* : A file describing the genes in the discovery network (can be loaded with gephi or cytoscape)

 Format:

Id | query | label | iteration
-- | ----- | ----- | ---------
node1 | query3 | "QuerySet:182893" | 0
node2 | query3 | "93_Amamu1:182893" | 1
node3 | query3 | "92_Aspsy1:40389" | 1

 * *outdir/run/graph/edges.tsv* : A file describing the discovery order in the discovery network

### Running your own problem

Please look at the example JSON file, and the default parameters in (pipeline_components/defaults.json)
The general format takes this form:

```json
{
  "genomes" : {  },
  "queries" : {  },
  "outdir" : "./output",
}
```

You can modify a number of parameters.
These can be found in (pipeline_components/defaults.json).

You can validate your configuration file with `ggMatch -v config.json`

#### Adding genomes

A genome is a set of proteins.
Provide for each genome a multisequence fasta file of protein sequences.
For each genome, add the location of the fasta file in the JSON file, indexed by a <genomeID>, under the "genomes" heading:
For example:
``` json
"genomes": {
  "genome1" : { "prots" : "proteins_genome1.fasta" },
  "genome2" : { "prots" : "proteins_genome2.fasta" }
}
```

#### Defining a query

A query takes the form of a multisequence fasta file.
A simple query is a single seed sequence, but if you have pre-existing knowledge of another ortholog, you can provide multiple seed sequences in the same query file.

If the genomes from which these queries originate also exist in your genome, you can prefix the sequence description in the fasta file with "genomeID:".
This will link the query sequence to the genome.
If that query is identified as a match against any other genome, then the reciprocal blast will be performed against that genome, rather than the set of query sequences.

For example:

```
>genome1:gene1
MPDDVWSGSSTCSLSSDGMSVRKDMKPEFHRAWPRCTAKAMDLEINEKMPHNETTEVAGVTKIKAVEAVG
GKTGKYIMYAGLAMVMVIYELDNSTVGTYRNFASSDFHQLGKLATLNTAASIITAIFKPPIAKLSDVLGR
GEAYVVTLTFYILSYILC
>prot2
MVAHNFSPRDAQFLTYTNGVSQALMGMGTGLLMYRYRTYKWIGVAGAVIRLVGYGVMVRLRTNESSIAEL
FIVQLVQGIGSGIIETIIIVAAQISVPHAELAQVTSLVMLGTFLGNGIGSAVAGAIYTNQLRDRLEIHLG
PGAAEGQLATLYNSITDRLPEWGTAERTAVNQALGDGHNLVQVTPDSSRSDSLDIEKPKARCF
```

These queries need to be added to the configuration JSON file:

```
  "queries" : {
    "query1" : "query_prots.fasta"
  }
```


