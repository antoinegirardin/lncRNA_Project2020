#!/usr/bin/env bash


#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=03:00:00
#SBATCH --job-name=fastqc
#SBATCH --mail-user=antoine.girardin@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/agirardin/output/output_fastqc_%j.o
#SBATCH --error=/data/users/agirardin/error/error_fastqc_%j.e

module add UHTS/Quality_control/fastqc/0.11.7
module add UHTS/Analysis/MultiQC/1.8;

READS_DIR=/data/courses/rnaseq/lncRNAs/Project2
HOME_DIR=/data/users/agirardin
mkdir -p $HOME_DIR/fastQC

cd $READS_DIR/fastq

FILES=(*.fastq.gz)

cd $HOME_DIR/fastQC

for i in `seq 0 23`
do fastqc -t 6 $READS_DIR/fastq/${FILES[${i}]} -o $HOME_DIR/fastQC
done 

multiqc .