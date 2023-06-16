#!/bin/bash
#SBATCH --mem 32G
#SBATCH -c 8
#SBATCH -t 14-00:00

set -eou pipefail

reads_dir=$1 # Where the gzipped fastq files are (basename of directory is used as sample name)
ref_fasta=$2 # Reference genome to align to (must be indexed before running this script)
gff=$3       # Annotations for reference genome in GFF format
chromosome_lengths=$4 # Path to csv file of chromosome lengths
out_dir=$5   # Directory in which to create a results directory with the name of the sample

source ~/mambaforge/bin/activate envs/pgt-pgs-markers

sample=$(basename $reads_dir)
out_dir=$out_dir/$sample

mkdir -p $out_dir/trimmed

for fastq in $reads_dir/*.fastq.gz; do
    gunzip -c $fastq | fastx_trimmer -f14 -Q33 >$out_dir/trimmed/$(basename ${fastq%.gz}) &
done
wait

prefix=$out_dir/$sample

bwa-0.7.5/bwa mem -M -t 8 $ref_fasta $out_dir/trimmed/* | samtools view -S -b - >$prefix.bam
samtools sort -@ 8 $prefix.bam >$out_dir/"$sample"_sorted.bam
samtools index $out_dir/"$sample"_sorted.bam
samtools mpileup -f $ref_fasta $out_dir/"$sample"_sorted.bam > $prefix.pileup
gzip $out_dir/trimmed/* &

python3 src/pileup_to_coverage_tables.py \
    --pileup $prefix.pileup --chromosome_lengths $chromosome_lengths --out_dir $out_dir &

python3 src/pileup_to_consensus.py \
    --pileup $prefix.pileup --ref $ref_fasta --out $prefix.fna
# Extract the coding sequence from the consensus
gffread -x $prefix.ffn -g $prefix.fna $gff 
python3 src/extract_every_third_base.py $prefix.ffn > $out_dir/"$sample"_third_bases.ffn
# Create a FASTA file with a single record containing every 3rd position in every codon
echo \>$sample > "$prefix"_third_bases_concat.ffn
cat "$prefix"_third_bases.ffn | grep -v \> | tr -d '\n' >> "$prefix"_third_bases_concat.ffn
echo >> "$prefix"_third_bases_concat.ffn

gzip $prefix.pileup &

wait
exit 0
