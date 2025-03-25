#!/bin/bash
#$ -N knitR2html
#$ -S /bin/bash
#$ -l tmem=15G
#$ -l h_vmem=15G
#$ -l h_rt=01:00:00
#$ -wd /SAN/ghlab/pophistory/Alice/hvCpG_project/code/2024_hvCpG/logs # one err and out file per sample

# Usage: qsub knitR2html.sh <filename.R>

if [ -z "$1" ]; then
    echo "Error: No filename provided"
    echo "Usage: qsub $0 <filename.R>"
    exit 1
fi

FILENAME="$1"
[[ $FILENAME == *.R ]] || FILENAME="${FILENAME}.R"

if [ ! -f "$FILENAME" ]; then
    echo "Error: File $FILENAME not found"
    exit 1
fi

# Be sure to have pandoc running
export PATH=$HOME/bin:$PATH

# Directly execute Rscript command
Rscript -e "rmarkdown::render('$FILENAME', output_format = 'html_document')"
