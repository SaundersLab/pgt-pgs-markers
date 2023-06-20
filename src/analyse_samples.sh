#!/bin/bash
#SBATCH --mem 8G
#SBATCH -J analyse_samples
#SBATCH -c 1
#SBATCH -o logs/analyse_samples.out
#SBATCH -e logs/analyse_samples.err
#SBATCH -t 14-00:00

set -eou pipefail

partitions=$1

ref_fasta=data/reference/Pgt_201_B1_AssemblyScaffolds2.fasta
gff=data/reference/Pgt_201_B1_GeneCatalog_20191015.gff3
chromosome_lengths=data/reference/chromosome_lengths.csv

source ~/mambaforge/bin/activate envs/pgt-pgs-markers

mkdir -p logs/samples

python3 src/chromosome_lengths.py $ref_fasta $chromosome_lengths

bwa-0.7.5/bwa index $ref_fasta

dependency=afterok
for sample_reads_dir in data/reads/*; do
    sample=$(basename $sample_reads_dir)
    sample_reads_dir=data/reads/$sample
    dependency=$dependency:$(
        sbatch \
            --parsable -J $sample -p $partitions \
            -o logs/samples/$sample.out -e logs/samples/$sample.err \
            src/analyse_sample.sh $sample_reads_dir $ref_fasta $gff $chromosome_lengths results/samples
    )
done
sbatch --kill-on-invalid-dep yes -o /dev/null --wait -d $dependency --wrap ""

# Combine the coverage data from different samples, then calculate statistics
# that are used to define marker regions
mkdir -p results/coverage results/coverage_stats logs/coverage
dependency=afterok
for chromosome in chr_{1..18}; do
    dependency=$dependency:$(
        sbatch \
            --parsable -J $chromosome -p $partitions --mem 32G -c 1 \
            -o logs/coverage/$chromosome.out -e logs/coverage/$chromosome.err \
            --wrap "
                source ~/mambaforge/bin/activate envs/pgt-pgs-markers
                paste -d, results/samples/*/coverage/$chromosome.csv > results/coverage/$chromosome.csv
                python3 src/coverage_to_coverage_stats.py \
                    results/coverage/$chromosome.csv \
                    results/coverage_stats/$chromosome.csv
                gzip results/coverage/$chromosome.csv
            "
    )
done
sbatch --kill-on-invalid-dep yes -o /dev/null --wait -d $dependency --wrap ""
python3 src/coverage_stats_to_markers.py $ref_fasta results/coverage_stats results/markers.csv

mkdir -p results/tree
cat results/samples/*/*_third_bases_concat.ffn >results/tree/msa.ffn

sbatch -J tree -p $partitions --mem 400G -J raxml-ng -c 32 -o logs/raxml-ng.out -e logs/raxml-ng.err -t 14-00:00 \
    --wrap "
    source ~/mambaforge/bin/activate envs/pgt-pgs-markers
    raxml-ng \
        --all --msa results/tree/msa.ffn --model GTR+G --tree pars{10},rand{10} \
        --prefix results/tree/tree --bs-trees 100 --threads 32 --seed 9 --redo
    "

python3 src/tree_input_coverage_filter.py results/tree/msa.ffn 80 results/tree/msa_80pct_covered.ffn

jobid=$(
    sbatch --parsable -J raxml-ng_80pct_covered -p $partitions --mem 400G -c 32 -o logs/raxml-ng_80pct_covered.out \
    -e raxml-ng_80pct_covered.err -t 14-00:00 --wrap "
source ~/mambaforge/bin/activate envs/pgt-pgs-markers
raxml-ng \
    --all --msa results/tree/msa_80pct_covered.ffn --model GTR+G --tree pars{10},rand{10} \
    --prefix results/tree/msa_80pct_covered --bs-trees 100 --threads 32 --seed 9 --redo
"
)
sbatch --kill-on-invalid-dep yes -o /dev/null --wait -d afterok:$jobid --wrap ""
sbatch -p $partitions src/report.sh

mkdir -p figures
python3 src/draw_tree.py data/metadata/metadata.csv \
    results/tree/msa_80pct_covered.raxml.support \
    figures/tree
