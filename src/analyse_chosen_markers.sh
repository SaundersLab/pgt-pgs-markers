#!/bin/bash

results_dir=$1

mkdir -p figures

mkdir -p data/reference/rust_assembly_sequences
for assembly in GCA_002762355.2 GCA_008522505.1 GCF_000149925.1 GCA_008520325.1 GCA_903797515.1 ; do
    mkdir -p temp_extract_assemblies
    echo $assembly
    echo "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$assembly/download?include_annotation_type=GENOME_FASTA&filename=$assembly.zip"
    curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$assembly/download?include_annotation_type=GENOME_FASTA&filename=$assembly.zip" -H "Accept: application/zip"
    mv $assembly.zip temp_extract_assemblies/
    pushd temp_extract_assemblies
    unzip $assembly.zip
    popd
    mv temp_extract_assemblies/ncbi_dataset/data/$assembly/*.fna data/reference/rust_assembly_sequences/
    sleep 2
    rm -r temp_extract_assemblies
done

pushd data/reference/rust_assembly_sequences/
for assembly in $(ls *_*_*) ; do
    short_name=$(echo $assembly | cut -f1,2 -d_)
    mv $assembly $short_name.fa
done
popd

mkdir -p results/samples_msas
python src/draw_products_msas.py &
tail -n +2 data/primers/marker_primer_regions.csv | while IFS=, read -r chromosome start end ; do
    range=$chromosome:$start-$end
    consensus=results/samples_msas/"$chromosome"_samples_consensuses.fasta
    msa=results/samples_msas/"$chromosome"_msa.fasta
    > $consensus 
    for bam in $(ls $results_dir/samples/*/*_sorted.bam) ; do
        sample=$(basename $bam | cut -d_ -f1)
        [ $sample == SAMEA104219957 ] || [ $sample == SAMEA104219958 ] && continue
        echo ">$sample" >> $consensus        
        samtools consensus -d 2 -l 10000 -A -r $range $bam | tail -n +2 >> $consensus
    done
    muscle -align $consensus -output $msa
done
python src/draw_sample_msas.py
wait < <(jobs -p)
