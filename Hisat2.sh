#! /usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=30:00:00
#SBATCH --job-name=hisat2
#SBATCH --mail-user=antoine.girardin@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/agirardin/output/output_hisat2_%j.o
#SBATCH --error=/data/users/agirardin/error/error_hisat2_%j.e

module add UHTS/Aligner/hisat/2.2.1
module add UHTS/Analysis/samtools/1.10;

READS_DIR=/data/courses/rnaseq/lncRNAs/Project2/
HOME_DIR=/data/users/agirardin
mkdir -p $HOME_DIR/hisat2

## Get reference genome indexes from https://daehwankimlab.github.io/hisat2/download/#h-sapiens
# wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# tar -xzf grch38_genome.tar.gz

cd $READS_DIR/fastq
FILES=(*.fastq.gz)

cd $HOME_DIR/hisat2
ln -s $READS_DIR/reference-files/grch38 .

for i in `seq 0 11`
do hisat2 -x grch38/genome -1 $READS_DIR/fastq/${FILES[2*${i}]} -2 $READS_DIR/fastq/${FILES[2*${i}+1]} -S Replicate_${i}.sam
samtools view -bS Replicate_${i}.sam > Replicate_${i}.bam
done