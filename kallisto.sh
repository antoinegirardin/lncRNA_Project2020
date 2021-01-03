#!/usr/bin/env bash


#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12000M
#SBATCH --time=24:00:00
#SBATCH --job-name=kallisto
#SBATCH --mail-user=antoine.girardin@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/agirardin/output/output_kallisto_%j.o
#SBATCH --error=/data/users/agirardin/error/error_kallisto_%j.e

## Load Kallisto and cufflinks modules
module add UHTS/Analysis/kallisto/0.46.0;
module add UHTS/Assembler/cufflinks/2.2.1;

## Define working directories
READS_DIR=/data/courses/rnaseq/lncRNAs/Project2
HOME_DIR=/data/users/agirardin
mkdir -p $HOME_DIR/kallisto

## Move to reads directory and create files list
cd $READS_DIR/fastq/
FILES=(*.fastq.gz)

## Move to results directory
cd $HOME_DIR/kallisto

## Read reference files
## GRCh38.p13.genome.fa : ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh38.p13.genome.fa.gz
## meta-assembly.gtf created at with StringTie.sh
ln -s $READS_DIR/reference-files/GRCh38.p13.genome.fa GRCh38.fa 
ln -s $HOME_DIR/meta-assembly.gtf .

## Create .fasta and index from transcriptome assembly and reference genome
gffread -w meta-assembly.fa -g GRCh38.fa meta-assembly.gtf
kallisto index -i meta-assembly.idx meta-assembly.fa

## Loop over each replicates and quantify the expression with kallisto quant
for i in `seq 0 11`
do kallisto quant -i meta-assembly.idx -b 100 -o $HOME_DIR/kallisto/results/Replicate_${i} $READS_DIR/fastq/${FILES[2*${i}]} $READS_DIR/fastq/${FILES[2*${i}+1]}
done


