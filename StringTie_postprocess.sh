#! /usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=00:10:00
#SBATCH --job-name=StringTie
#SBATCH --mail-user=antoine.girardin@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/agirardin/output/output_StringTie_%j.o
#SBATCH --error=/data/users/agirardin/error/error_StringTie_%j.e

cd StringTie

awk '{if($3=="transcript"){print $0}}' meta-assembly.gtf > meta-assembly_transcripts.gtf
awk '{if($10~ /MSTRG/){print $0}}' meta-assembly_transcripts.gtf > meta-assembly_noveltranscripts.gtf
awk '{if($13~ /gene/){print $0}}' meta-assembly_noveltranscripts.gtf > meta-assembly_notannotatedtranscripts.gtf
awk '{if($13 !~ /gene/){print $0}}' meta-assembly_transcripts.gtf > meta-assembly_noveltranscripts.gtf
awk '{if($10~ /ENSG/){print $0}}' meta-assembly_transcripts.gtf > meta-assembly_annotatedtranscripts.gtf
awk '{print $14}' meta-assembly_annotatedtranscripts.gtf | sort | uniq -c > meta-assembly_annotatedgenes.txt

wc -l meta-assembly_transcripts.gtf > meta-assembly_summary.csv
wc -l meta-assembly_annotatedtranscripts.gtf >> meta-assembly_summary.csv
wc -l meta-assembly_noveltranscripts.gtf >> meta-assembly_summary.csv
wc -l meta-assembly_notannotatedtranscripts.gtf >> meta-assembly_summary.csv
wc -l meta-assembly_annotatedgenes.txt >> meta-assembly_summary.csv

awk 'BEGIN{k=""}{if($3=="transcript"){print i,j,k; i=0; j=$5-$4; k=$12}if($3=="exon"){i+=1}}END{print i,j,k}' meta-assembly.gtf > meta-assembly_count.csv
sed -i '1d' meta-assembly_count.csv
sed -i -e 's/\"//g' meta-assembly_count.csv
sed -i -e 's/;//g' meta-assembly_count.csv
awk 'BEGIN{i=j=k=0}{i+=1;j+=$1;if($1==1){k+=1}}END{print i,",",j,",",k}' meta-assembly_count.csv > meta-assembly_summary1.csv
