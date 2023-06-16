#!/bin/bash

set -eou pipefail

reads_dir=$1

tail -n +2 data/metadata/reads.csv | while IFS=, read -r run_accession sample_accession sample_name data_source ebi_fastq_ftp_r1 ebi_fastq_ftp_r2; do
    mkdir -p $reads_dir/$sample_accession
    echo "$run_accession"
    if [[ $data_source = "NCBI SRA" ]]; then
        fasterq-dump -f --progress --temp $reads_dir --outdir $reads_dir/$sample_accession --split-files $run_accession
    elif [[ $data_source = "EMBL-EBI ENA" ]]; then
        curl -O --output-dir $reads_dir/$sample_accession $ebi_fastq_ftp_r1
        curl -O --output-dir $reads_dir/$sample_accession $ebi_fastq_ftp_r2
    fi
done
gzip $reads_dir/*/*.fastq
