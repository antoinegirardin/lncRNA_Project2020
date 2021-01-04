# lncRNA_Project2020

Created by Antoine Girardin on december 2020.

This repository contains all scripts related to data treatment for RNA-seq course's project. In the folders you can find some results files that are not too heavy. For *.bam* and *.gtf* files or kallisto results, please go on the server: 

- binfservms01.unibe.ch:/data/users/agirardin

## Step 1: Quality control

The first step of this project is to control the quality of the reads with *fastQC.sh* script. The *.html* results can be found in the folder */FastQC*. Please check multiqc_report.html for integrated analysis.

## Step 2: Read mapping

*Hisat2.sh* is used for mapping the read with a reference genome. It creates one *.sam* file per replicates. It also converts the *.sam* file into a *.bam* file. These files are not available on github but you can check the summary files in */Hisat2* folder.

To pursue to the next step, the script *sortSam.sh* is used to sort and merge replicates' *.bam* files of the same cell subtype. These are called as follow:

- Holoclone:    Sample1.bam
- Meroclone:    Sample2.bam
- Paraclone:    Sample3.bam
- Parental:     SampleP.bam

## Step 3: Transcriptome assembly

For transcriptome assembly, *StringTie.sh* is used and it create a *.gtf* file for each *.bam* file created at the previous step. The 4 samples files are merge in one *meta-assembly.gtf* file. 

The *meta-assembly.gtf* is processed with *awk* in *StringTie_postprocess.sh* to count the number of each type of transcripts. It creates new filtered version of *meta-assembly.gtf*. All *.gtf* files are on the server but a summary is located in */StringTie* folder.

*Cpat.sh* and *intersect.sh* are used to determine the quality of each transcript from *meta-assembly-transcript.gtf* file. Results can be found respectively in */Cpat* and */intersect* folders.

## Step 4: Quantification

The quantification is done with *kallisto.sh* that generates a results folder per replicates. A summary of these results is available in */kallisto* folder. This summary is produced with *kallisto_postprocess.sh* script.

## Step 5: Differential expression

Differential expression is calculated with *Sleuth.R* that uses *kallisto*'s results and other intermediate results such as:

- cpat_output
- cage_peak_intersect.txt
- polyAsite_intersect.txt
- Intergenic_intersect.txt
- meta-assembly_count.csv

## Report

The results are presented in the report named *RNA-seq_Report_ANG.html*.

For further questions, please contact me at the following e-mail address: antoine.girardin@unifr.ch