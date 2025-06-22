#WD=/media/inter/mkapun/projects/MuseomicsWorkshop2025
WD=$1

eval "$(conda shell.bash hook)"

conda create \
    -p ${WD}/scripts/programs \
    -y \
    -c bioconda \
    -c conda-forge \
    mamba

conda activate ${WD}/scripts/programs

## (1) install requirements for the Data Preparation scripts
# mamba create \
#     -p ${WD}/scripts/prepareData \
#     -y \
#     -c bioconda \
#     -c conda-forge \
#     fastp seqkit=2.10 bwa samtools pigz

## (2a) install requirements for the Workshop scripts

mamba install \
    -y \
    -c bioconda \
    -c conda-forge \
    fastp seqkit minimap2 bwa samtools blast bcftools r r-base r-tidyverse bioconductor-ggtree r-phangorn r-phytools r-ape r-patchwork

conda deactivate

## (2b) install mapdamage2 in a separate environment
mamba create \
    -p ${WD}/scripts/mapdamage2 \
    -y \
    -c bioconda mapdamage2
