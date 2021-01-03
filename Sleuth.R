#!/usr/bin/env Rscript

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10000M
#SBATCH --time=03:00:00
#SBATCH --job-name=sleuth
#SBATCH --mail-user=antoine.girardin@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/agirardin/output/output_sleuth_%j.o
#SBATCH --error=/data/users/agirardin/error/error_sleuth_%j.e
#SBATCH --array=0

# module add R/3.6.1;

setwd("/data/users/agirardin/sleuth")
rm(list=ls())

suppressPackageStartupMessages(library("sleuth"))
suppressPackageStartupMessages(library("biomaRt"))

sample_id <- dir(file.path("..","kallisto","results"))
kal_dirs <- file.path("..","kallisto","results", sample_id)
s2c <- read.table(file.path("exp_design.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

cat("s2c \n")

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

cat("t2g \n")

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE,target_mapping = t2g)

cat("so \n")

so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

cat("results \n")

## fold change
sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm') # 'obs_raw', 'est_counts'
## /!\ Add 0.1 to avoid small value rounded to zero!!!
sleuth_matrix <- sleuth_matrix + 0.1
sleuth_table <- cbind(sleuth_table,fold_change1=rep(0,dim(sleuth_table)[1]))
sleuth_table <- cbind(sleuth_table,fold_change2=rep(0,dim(sleuth_table)[1]))
sleuth_table <- cbind(sleuth_table,fold_change3=rep(0,dim(sleuth_table)[1]))
sleuth_table <- cbind(sleuth_table,logqvalue=-log10(sleuth_table$qval))
fold_change <- data.frame(geneName=rownames(sleuth_matrix),fold_change1=log2(rowMeans(sleuth_matrix[,1:3])/rowMeans(sleuth_matrix[,10:12])),fold_change2=log2(rowMeans(sleuth_matrix[,4:6])/rowMeans(sleuth_matrix[,10:12])),fold_change3=log2(rowMeans(sleuth_matrix[,7:9])/rowMeans(sleuth_matrix[,10:12])))
cat("fold \n")

## coding potential
cpat_score <- read.table("/data/users/agirardin/cpat/cpat_output")
cpat_score <- cbind(cpat_score, name=rownames(cpat_score))
sleuth_table <- cbind(sleuth_table,coding_prob=rep(2,dim(sleuth_table)[1]))
cat("cpat \n")

## transcript quality
metaAssemblyCount <- read.table("/data/users/agirardin/StringTie/meta-assembly_count.csv",header = F)
sleuth_table <- cbind(sleuth_table,nb_exon=rep(0,dim(sleuth_table)[1]))
sleuth_table <- cbind(sleuth_table,base_length=rep(0,dim(sleuth_table)[1]))
cat("exon \n")

## 5' end quality
cage_peak_intersect <- read.table("/data/users/agirardin/intersect/cage_peak_intersect.txt",header = F)
sleuth_table <- cbind(sleuth_table,cage_peak=rep(0,dim(sleuth_table)[1]))
cat("5 \n")

## 3' quality
polyAsite_intersect <- read.table("/data/users/agirardin/intersect/polyAsite_intersect.txt",header = F)
sleuth_table <- cbind(sleuth_table,polyAsite=rep(0,dim(sleuth_table)[1]))
cat("3 \n")

## Intergenic
Intergenic_intersect <- read.table("/data/users/agirardin/intersect/Intergenic_intersect.txt",header = F)
sleuth_table <- cbind(sleuth_table,intergenic=rep(0,dim(sleuth_table)[1]))
cat("intergenic \n")

for (i in 1:dim(sleuth_table)[1]){
  sleuth_table[i,15:17] <- fold_change[fold_change$geneName==sleuth_table$target_id[i],2:4]
  sleuth_table$coding_prob[i] <- cpat_score$coding_prob[cpat_score$name==sleuth_table$target_id[i]]
  sleuth_table[i,20:21] <- metaAssemblyCount[metaAssemblyCount[,3]==sleuth_table$target_id[i],1:2]
  sleuth_table$cage_peak[i] <- sum(cage_peak_intersect[,7]==sleuth_table$target_id[i])
  sleuth_table$polyAsite[i] <- sum(polyAsite_intersect[,7]==sleuth_table$target_id[i])
  sleuth_table$intergenic[i] <- sum(Intergenic_intersect[,7]==sleuth_table$target_id[i])
}

cat("integrative analysis \n")

sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.1)
sleuth_significant_codingprot <- dplyr::filter(sleuth_significant, coding_prob <= 0.364)
sleuth_significant_foldchange1 <- dplyr::filter(sleuth_significant_codingprot, abs(fold_change1) >= 1)
sleuth_significant_foldchange2 <- dplyr::filter(sleuth_significant_codingprot, abs(fold_change2) >= 1)
sleuth_significant_foldchange3 <- dplyr::filter(sleuth_significant_codingprot, abs(fold_change3) >= 1)
# head(sleuth_significant_codingprot, 20)

save.image()

# sleuth_live(so)




