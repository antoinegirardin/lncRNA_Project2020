#! /usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20000M
#SBATCH --time=01:00:00
#SBATCH --job-name=StringTie
#SBATCH --mail-user=antoine.girardin@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/agirardin/output/output_StringTie_%j.o
#SBATCH --error=/data/users/agirardin/error/error_StringTie_%j.e

## Load StringTie module
module add UHTS/Aligner/stringtie/1.3.3b;

## Define working directory and move to results directory
READS_DIR=/data/courses/rnaseq/lncRNAs/Project2/
HOME_DIR=/data/users/agirardin
mkdir -p $HOME_DIR/StringTie
cd $HOME_DIR/StringTie

## Load .bam file 
ln -s $HOME_DIR/hisat2/Sample*.bam .

## Load reference annotation file from link below
## wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.chr_patch_hapl_scaff.annotation.gtf.gz
## Reference file strored in the following location
ln -s $READS_DIR/reference-files/gencode.v35.chr_patch_hapl_scaff.annotation.gtf .

## Launch StringTie on each sample and merge all .gtf file to a meta-assembly.gtf file
stringtie Sample1.bam -G gencode.v35.chr_patch_hapl_scaff.annotation.gtf -o Sample1.gtf
stringtie Sample2.bam -G gencode.v35.chr_patch_hapl_scaff.annotation.gtf -o Sample2.gtf
stringtie Sample3.bam -G gencode.v35.chr_patch_hapl_scaff.annotation.gtf -o Sample3.gtf
stringtie SampleP.bam -G gencode.v35.chr_patch_hapl_scaff.annotation.gtf -o SampleP.gtf
stringtie --merge -G gencode.v35.chr_patch_hapl_scaff.annotation.gtf Sample1.gtf Sample2.gtf Sample3.gtf SampleP.gtf -o meta-assembly.gtf
