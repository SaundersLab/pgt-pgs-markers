#!/bin/bash
#SBATCH --mem 120G
#SBATCH -c 20
#SBATCH -t 14-00:00
#SBATCH -o logs/report.out
#SBATCH -e logs/report.err

set -eou pipefail

source ~/mambaforge/bin/activate envs/pgt-pgs-markers

mkdir -p reports/alignment

for bam in results/samples/*/*_sorted.bam ; do
    sample=$(basename $bam | cut -d_ -f1)
    echo $sample
    samtools flagstat -O json $bam > reports/alignment/$sample.json &
done
wait

python src/report_table.py

