# Reproducing the Analysis

Analysis was run in 3 different environments:

- HPC "compute nodes" with no internet access
- HPC "software node" with limited internet access
- laptop running MacOS with full internet access and access to storage space visible to HPC nodes

## Setup

The following was run on a "software node":

```bash
# mamba forge was installed into home directory and sourced
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
source ~/mambaforge/bin/activate

# environment was created from the pgt-pgs-markers/envs directory
git clone https://github.com/SaundersLab/pgt-pgs-markers.git
cd pgt-pgs-markers/envs
mamba env create -f pgt-pgs-markers-env.yml --prefix pgt-pgs-markers
cd ..

# installation of bwa 0.7.5 in the pgt-pgs-markers directory
curl -L -O https://github.com/lh3/bwa/archive/refs/tags/0.7.5.tar.gz
tar -xvf 0.7.5.tar.gz
cd bwa-0.7.5
make
```

The following was run on a laptop:

```bash
# mamba forge was installed into home directory and sourced
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
source ~/mambaforge/bin/activate

# environment was created from the pgt-pgs-markers/envs directory
git clone https://github.com/SaundersLab/pgt-pgs-markers.git
cd pgt-pgs-markers/envs
mamba env create -f download-reads-env.yml --prefix download-reads-env
conda activate $PWD/download-reads-env
cd ..

path_to_hpc_storage_from_laptop=/path/to/hpc/storage/from/laptop

# script to download reads was run
./src/download_reads.sh "$path_to_hpc_storage_from_laptop"/pgt-pgs-markers/data/reads
```

The reference genome and annotations were downloaded from <https://genome.jgi.doe.gov/portal/Pgt_201_B1/Pgt_201_B1.download.html> (account required)

- Pgt_201_B1_GeneCatalog_20191015.gff3.gz
- Pgt_201_B1_AssemblyScaffolds2.fasta.gz

These were unzipped and copied into `pgt-pgs-markers/data/reference`

The empty lines (882866 and 1018814) were removed from Pgt_201_B1_AssemblyScaffolds2.fasta.

## Alignments, Marker Finding, and Tree Building

The `partitions` variable was set to a comma separated list of partition names specific to the cluster.

```bash
mkdir -p logs
partitions=comma,separated,list,of,partition,names
sbatch -p $partitions src/analyse_samples.sh $partitions
```

## Primer design

Results from `results/markers.csv` were used with [Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/).

The following changes were made from the default parameters:

|Parameter              |Default|Value Used|
|-----------------------|-------|----------|
|# of primers to return |10     |30        |
|Primer Size Max        |25     |30        |
|GC clamp               |0      |1         |

The following assemblies were used (links below are to the Primer-BLAST page for each assembly):

*Pgt* assemblies:

- [21-0](https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?PRIMER_BLAST_SPEC=Assembly&PRIMER_SPECIFICITY_DATABASE_UID=105670684)
- [CRL 75-36-700-3](https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?PRIMER_BLAST_SPEC=Assembly&PRIMER_SPECIFICITY_DATABASE_UID=17203964)
- [Ug99](https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?PRIMER_BLAST_SPEC=Assembly&PRIMER_SPECIFICITY_DATABASE_UID=105649704)
- [UK-01](https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?PRIMER_BLAST_SPEC=Assembly&PRIMER_SPECIFICITY_DATABASE_UID=131661594)
- [99KS76A-1](https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?PRIMER_BLAST_SPEC=Assembly&PRIMER_SPECIFICITY_DATABASE_UID=66690264)

Wheat:

- [Chinese Spring](https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?PRIMER_BLAST_SPEC=Assembly&PRIMER_SPECIFICITY_DATABASE_UID=79283754)

Rye:

- [Lo7](https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?PRIMER_BLAST_SPEC=Assembly&PRIMER_SPECIFICITY_DATABASE_UID=31433893)

Primers were designed on the 21-0 assembly and then checked against all the other assemblies listed.

## Analysis of chosen markers

```bash
# assuming we're still in the pgt-pgs-markers directory and 
# download-reads-env environment is still active
conda deactivate
cd envs
mamba env create -f analyse-chosen-markers.yml --prefix analyse-chosen-markers
conda activate $PWD/analyse-chosen-markers
cd ..

./src/analyse_chosen_markers.sh "$path_to_hpc_storage_from_laptop"/pgt-pgs-markers/results
```

Figures produced were then rearranged in a vector graphics editor.
