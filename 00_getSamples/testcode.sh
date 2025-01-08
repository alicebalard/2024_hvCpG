#$ -l tmem=2G
#$ -l h_vmem=2G
#$ -l h_rt=1:0:0
#$ -S /bin/bash
#$ -j y
#$ -N MyTESTJOBNAME

#/share/apps/genomics/bowtie2-2.4.1/bowtie2 mysamples genome -n 20

echo "Hello Marina"
