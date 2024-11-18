## DL files from a SRA list
#/share/apps/genomics/sratoolkit.3.0.2/bin/prefetch --max-size 30GB --option-file sralist

## Dataset1
/share/apps/genomics/sratoolkit.3.0.2/bin/prefetch --max-size 30GB --option-file sraAccList1.txt


/share/apps/genomics/sratoolkit.3.0.2/bin/fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ~/ncbi/public/sra/SRR2121685.sra



