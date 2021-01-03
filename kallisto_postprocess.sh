#!/usr/bin/env bash


#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00:02:00
#SBATCH --job-name=kallisto
#SBATCH --mail-user=antoine.girardin@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/agirardin/output/output_kallisto_%j.o
#SBATCH --error=/data/users/agirardin/error/error_kallisto_%j.e

cd kallisto/results

for i in `ls`
do echo $i >> ../kallisto_summary.txt
awk 'BEGIN{i=0}{if($5>0){i+=$5}}END{print i}' ${i}/abundance.tsv >> ../kallisto_summary.txt
awk 'BEGIN{i=j=0}{i+=1;if($5>0){j+=1}}END{k=100*j/i;print k}' ${i}/abundance.tsv >> ../kallisto_summary.txt
awk '$1~/p_pseudoaligned/{print $2}' ${i}/run_info.json >> ../kallisto_summary.txt
done

cd ..
awk 'BEGIN{i=0;j=""}{i+=1;j=j" "$0;if(i==4){i=0;print j; j=""}}' kallisto_summary.txt > kallisto_summary1.txt
sed -i -e 's/,//g' kallisto_summary1.txt