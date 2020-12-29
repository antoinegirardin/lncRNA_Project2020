#!/usr/bin/env bash


#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:05:00
#SBATCH --job-name=Intersect
#SBATCH --mail-user=antoine.girardin@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/agirardin/output/output_intersect_%j.o
#SBATCH --error=/data/users/agirardin/error/error_intersect_%j.e

module add UHTS/Analysis/BEDTools/2.29.2;
module add Utility/UCSC-utils/359;

READS_DIR=/data/courses/rnaseq/lncRNAs/Project2
HOME_DIR=/data/users/agirardin
mkdir -p $HOME_DIR/intersect
cd $HOME_DIR/intersect

ln -s $READS_DIR/reference-files/Integrative_analysis/hg19.cage_peak_phase1and2combined_coord.bed cage_peak.bed
ln -s $READS_DIR/reference-files/Integrative_analysis/hg19ToHg38.over.chain.gz .
ln -s $READS_DIR/reference-files/Integrative_analysis/atlas.clusters.2.0.GRCh38.96.bed polyAsite.bed
ln -s $READS_DIR/reference-files/gencode.v35.chr_patch_hapl_scaff.annotation.gtf .
ln -s $HOME_DIR/meta-assembly.gtf .

liftOver cage_peak.bed hg19ToHg38.over.chain.gz cage_peak_hg38.bed unMapped

bedtools intersect -wa -a meta-assembly.gtf -b cage_peak_hg38.bed > cage_peak_intersect.txt
bedtools intersect -wa -a meta-assembly.gtf -b polyAsite.bed > polyAsite_intersect.txt
bedtools intersect -v -a meta-assembly.gtf -b gencode.v35.chr_patch_hapl_scaff.annotation.gtf > Intergenic_intersect.txt

