#!/usr/bin/env Rscript

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10000M
#SBATCH --time=00:30:00
#SBATCH --job-name=sleuth
#SBATCH --mail-user=antoine.girardin@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/agirardin/output/output_sleuth_%j.o
#SBATCH --error=/data/users/agirardin/error/error_sleuth_%j.e
#SBATCH --array=0

# module add R/3.6.1;

setwd("/data/users/agirardin/sleuth")

suppressPackageStartupMessages(library("sleuth"))
suppressPackageStartupMessages(library("biomaRt"))

sample_id <- dir(file.path("..","kallisto","results"))
kal_dirs <- file.path("..","kallisto","results", sample_id)
s2c <- read.table(file.path("exp_design.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE,target_mapping = t2g)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
# models(so)
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
# head(sleuth_significant, 20)

sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm') # 'obs_raw', 'est_counts'
## /!\ Add 0.1 to avoid small value rounded to zero!!!

sleuth_table <- cbind(sleuth_table,fold_change1=rep(0,dim(sleuth_table)[1]))
sleuth_table <- cbind(sleuth_table,fold_change2=rep(0,dim(sleuth_table)[1]))
sleuth_table <- cbind(sleuth_table,fold_change3=rep(0,dim(sleuth_table)[1]))

fold_change <- data.frame(geneName=rownames(sleuth_matrix),fold_change1=log2(rowMeans(sleuth_matrix[,1:3])/rowMeans(sleuth_matrix[,10:12])),fold_change2=log2(rowMeans(sleuth_matrix[,4:6])/rowMeans(sleuth_matrix[,10:12])),fold_change3=log2(rowMeans(sleuth_matrix[,7:9])/rowMeans(sleuth_matrix[,10:12])))

for (i in 1:dim(sleuth_table)[1]){
  sleuth_table[i,15:17] <- fold_change[fold_change$geneName==sleuth_table$target_id[i],2:4]
}

sleuth_table <- cbind(sleuth_table,logqvalue=-log10(sleuth_table$qval))

save.image()

# sleuth_live(so)




