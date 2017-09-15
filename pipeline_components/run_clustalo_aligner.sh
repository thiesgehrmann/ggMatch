#!/usr/bin/env sh

SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

configfile="$1"
threads="$2"

snakemake --snakefile "$SCRIPTDIR/clustalo.Snakefile" --configfile "$configfile" --cores $threads --nolock --conda-prefix "$SCRIPTDIR/../.snakemake/conda" --use-conda
