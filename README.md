# gxsMatch
Greedy Cross-Species gene Matching tool

## Usage

### Dependencies

 - Snakemake
 - Conda

### Running the example

You can download and run the example dataset with the following commands:

```bash
  git clone https://github.com/thiesgehrmann/gxsMatch.git
  cd gxsMatch
  snakemake --nolock --configfile example/config.json --cores 20 --use-conda final
```

The --nolock argument is important, since I call a sub-instance of snakemake in the main snakemake file.

The example is defined by the 
