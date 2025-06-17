#!/bin/bash

###############################################################################
# Museomics Workshop 2025: Bioinformatics Pipeline
#
# This script describes the bioinformatics workflow for processing and analyzing
# sequencing data from the Museomics Workshop 2025.
#
# All explanatory text is included as commented blocks.
###############################################################################

################################################################################
# 1. Preparation
################################################################################

## System Requirements:
# - Linux/Unix (macOS) system with a BASH terminal
# - ~10GB free space
# - 16GB RAM / 4 cores
# - VScode installed: https://code.visualstudio.com/
# - Conda installed and in PATH: https://anaconda.org/anaconda/conda

# Set Working directory
WD="<your/path/to/MuseomicsWorkshop2025>"
# WD="/media/inter/mkapun/projects/MuseomicsWorkshop2025"

# Make sure conda and mamba are installed on your system
sh "${WD}/shell/requirements.sh" "${WD}"

################################################################################
# 2. Download and decompress the raw input FASTQ data
################################################################################

## Download and extract raw sequencing reads.
cd "${WD}/data" && mkdir -p "${WD}/data/raw_reads"

# Download raw reads
wget -O raw_reads.tar.gz "https://filesender.aco.net/download.php?token=18816b43-8209-46cf-ae04-40cece8dc2c1&files_ids=895410"

# Untar the compressed folder within ${WD}/data
tar -xzf raw_reads.tar.gz -C "${WD}/data/raw_reads"

## optionally download mapDamage2 results
mkdir ${WD}/results/mapDamage
wget -O mapDamage2.tar.gz "https://filesender.aco.net/download.php?token=cc26f12d-2afe-4126-8383-ff3597b5412c&files_ids=896942"
tar -xzf ${WD}/data/mapDamage2.tar.gz -C ${WD}/results/mapDamage

################################################################################
# 3. Datasets
################################################################################

: '
# Datasets table:

# | Library  | Name    | Age  | City     | Country | Wolbachia    | Type     | SRA        |
# |----------|---------|------|----------|---------|--------------|----------|------------|
# | DGRP370  | DGRP370 | 2003 | Raleigh  | USA     | wMel         | recent   | SRR834539  |
# | DGRP338  | DGRP338 | 2003 | Raleigh  | USA     | wMelCS       | recent   | SRR834513  |
# | ZI268    | ZI268   | 2003 | Ziawonga | Zambia  | wMel         | recent   | SRR189425  |
# | HG_15    | 19SL15  | 1933 | Lund     | Sweden  | unknown      | historic | SRR23876580|
# | HG0027   | 19SL3   | 1933 | Lund     | Sweden  | unknown      | historic | SRR23876574|
# | HG0029   | 18DZ5   | 1899 | Zealand  | Denmark | unknown      | historic | SRR23876565|
'

################################################################################
# 4. Preview Raw FASTQ Data
################################################################################

# Preview the first 10 lines of the forward read file for two datasets.
less "${WD}/data/raw_reads/ZI268_1.fastq.gz"
less "${WD}/data/raw_reads/18DZ5_1.fastq.gz"

# QUESTIONS:
# A) What is the structure of the datasets?
# B) How do the two files differ?

################################################################################
# 5. Trim Reads with fastp
################################################################################

# Trim and filter reads using fastp. Focus on a subset of 1,000,000 reads from each dataset.

mkdir -p "${WD}/data/trimmed_reads" && cd "${WD}/data/trimmed_reads"
conda activate "${WD}/scripts/programs"

# Loop through each library in datasets.csv and trim reads
while IFS="," read -r Library Name Age City Country Wolb Type SRA; do
    if [[ "${Library}" != "Library" ]]; then
        echo "Processing library ${Name}"
        fastp \
            -i "${WD}/data/raw_reads/${Name}_1.fastq.gz" \
            -I "${WD}/data/raw_reads/${Name}_2.fastq.gz" \
            -o "${WD}/data/trimmed_reads/${Name}_1_trimmed.fastq.gz" \
            -O "${WD}/data/trimmed_reads/${Name}_2_trimmed.fastq.gz" \
            --merge \
            --merged_out "${WD}/data/trimmed_reads/${Name}_merged.fastq.gz" \
            --length_required 25 \
            --dedup \
            --trim_poly_g \
            --html "${WD}/data/trimmed_reads/${Name}.html" \
            --json "${WD}/data/trimmed_reads/${Name}.json" \
            --detect_adapter_for_pe
    fi
done <"${WD}/data/datasets.csv"

# View HTML output files
firefox "${WD}/data/trimmed_reads/18DZ5.html"
firefox "${WD}/data/trimmed_reads/ZI268.html"

# QUESTIONS:
# A) What do the parameters mean?
# B) How do the read length distributions differ between the datasets?

################################################################################
# 6. Testing for Eukaryotic Contamination
################################################################################

# Download mitochondrial genomes for contamination screening and build a BLAST database.

mkdir -p "${WD}/data/refseq/contamination" && cd "${WD}/data/refseq/contamination"

# Download and rename mitochondrial genomes
wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=PP230540.1&db=nuccore&report=fasta&format=text" -O Gryllus_bimaculatus_mito.fasta
sed -i '' '1s/.*/>Gryllus_bimaculatus/' Gryllus_bimaculatus_mito.fasta

wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_012920.1&db=nuccore&report=fasta&format=text" -O Homo_sapiens_mito.fasta
sed -i '' '1s/.*/>Homo_sapiens/' Homo_sapiens_mito.fasta

wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_024511.2&db=nuccore&report=fasta&format=text" -O Drosophila_melanogaster_mito.fasta
sed -i '' '1s/.*/>Drosophila_melanogaster/' Drosophila_melanogaster_mito.fasta

wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=PV608434.1&db=nuccore&report=fasta&format=text" -O Anthrenus_verbasci_mito.fasta
sed -i '' '1s/.*/>Anthrenus_verbasci/' Anthrenus_verbasci_mito.fasta

# Concatenate all files into one
cat *_mito.fasta >mitochondrial.fasta

# Build nucleotide BLAST database
makeblastdb -in mitochondrial.fasta -dbtype nucl -out mitoDB

# Convert trimmed FASTQ files to FASTA and run BLAST against the mitochondrial database.

# Prepare results folder and output file
mkdir -p "${WD}/results/contamination"
>"${WD}/results/contamination/BLAST.tsv"

# Loop through all files in dataset, convert FASTQ to FASTA, then BLAST
while IFS="," read -r Library Name Age City Country Wolb Type SRA; do
    if [[ "${Library}" != "Library" ]]; then
        gunzip -c "${WD}/data/trimmed_reads/${Name}_1_trimmed.fastq.gz" |
            awk 'NR % 4 == 1 {print ">" substr($0, 2)} NR % 4 == 2 {print $0}' >"${WD}/data/trimmed_reads/${Name}_trimmed.fasta"

        gunzip -c "${WD}/data/trimmed_reads/${Name}_merged.fastq.gz" |
            awk 'NR % 4 == 1 {print ">" substr($0, 2)} NR % 4 == 2 {print $0}' >>"${WD}/data/trimmed_reads/${Name}_trimmed.fasta"

        blastn \
            -num_threads 4 \
            -evalue 1e-10 \
            -max_target_seqs 1 \
            -outfmt '6 qseqid sseqid slen qlen pident length mismatch evalue bitscore' \
            -db "${WD}/data/refseq/contamination/mitoDB" \
            -query "${WD}/data/trimmed_reads/${Name}_trimmed.fasta" |
            awk -v Name="${Name}" '$5>90 {print Name"\t"$0}' >>"${WD}/results/contamination/BLAST.tsv"
    fi
done <"${WD}/data/datasets.csv"

# View tabular output file
less "${WD}/results/contamination/BLAST.tsv"

# QUESTIONS:
# A) What do the columns mean?
# B) Do you see non-Drosophila hits?

# Plot the proportion of reads mapped to mitochondrial genomes in R
"${WD}/scripts/programs/bin/Rscript" -e "
library(tidyverse)
df <- read_tsv('${WD}/results/contamination/BLAST.tsv', col_names = FALSE)
df.prop <- df %>%
    select(1,3,5,9) %>%
    rename(Library = 1, Name = 2, ReadLength = 3, Evalue = 4) %>%
    group_by(Library,Name) %>%
    summarise(Count = n()) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    arrange(desc(Count))
plot <- ggplot(df.prop, aes(x = Library, y = Proportion, fill = Name)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(title = 'Proportion of Reads Mapped to Mitochondrial Genomes',
         x = 'Library',
         y = 'Proportion of Reads') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('${WD}/results/contamination/BLAST_proportion_plot.pdf', plot, width = 10, height = 6)
ggsave('${WD}/results/contamination/BLAST_proportion_plot.png', plot, width = 10, height = 6)
"

# QUESTIONS:
# A) What is the proportion of contamination?
# B) Are there differences between historic and recent datasets?

# Plot read length distributions for endogenous and exogenous DNA
"${WD}/scripts/programs/bin/Rscript" -e "
library(tidyverse)
df <- read_tsv('${WD}/results/contamination/BLAST.tsv', col_names = FALSE)
df.prop.RL <- df %>%
    select(1,3,5,9) %>%
    rename(Library = 1, Name = 2, ReadLength = 3, Evalue = 4) %>%
    group_by(Library,Name,ReadLength) %>%
    summarise(Count = n()) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    arrange(desc(Count))
plot.RL <- ggplot(df.prop.RL, 
    aes(x = ReadLength, 
        y = Count, 
        fill = Name)) +
    geom_bar(stat = 'identity', 
        position = 'dodge') +
    labs(title = 'Proportion of Reads Mapped to Mitochondrial Genomes by Read Length',
         x = 'Read Lengths',
         y = 'Read Counts') +
         facet_grid(Library~Name, scales = 'free_y') +
    guides(fill = guide_legend(title = 'Mitochondrial Genome')) +
    theme_bw() +
    scale_y_log10() +   
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) + theme(legend.position = 'bottom')
ggsave('${WD}/results/contamination/BLAST_proportion_plot_by_RL.pdf',
    plot.RL, 
    width = 7, 
    height = 5)
ggsave('${WD}/results/contamination/BLAST_proportion_plot_by_RL.png',
    plot.RL, 
    width = 7, 
    height = 5)
"

# QUESTIONS:
# A) How to interpret the differences in read lengths between the historic and recent samples?
# B) Do the read lengths differ between the exogenous and the endogenous DNA?
# C) What does this tell us about the contamination source?

################################################################################
# 7. Testing for DNA Degradation
################################################################################

# Map reads to reference genomes using minimap2 and analyze DNA damage patterns with mapDamage2.

mkdir -p "${WD}/results/minimap2" && cd "${WD}/results/minimap2"
while IFS="," read -r Library Name Age City Country Wolb Type SRA; do
    if [[ "${Library}" != "Library" ]]; then
        echo "Processing library ${Name}"
        minimap2 -ax sr --secondary=no -t 4 \
            "${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz" \
            "${WD}/data/raw_reads/${Name}_2L_100K_1_trimmed.fastq.gz" \
            "${WD}/data/raw_reads/${Name}_2L_100K_2_trimmed.fastq.gz" |
            samtools view -bS -F 4 - |
            samtools sort -o "${WD}/results/minimap2/${Name}_PE.bam"
        samtools index "${WD}/results/minimap2/${Name}_PE.bam"

        minimap2 -ax sr --secondary=no -t 4 \
            "${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz" \
            "${WD}/data/raw_reads/${Name}_2L_100K_merged.fastq.gz" |
            samtools view -bS -F 4 - |
            samtools sort -o "${WD}/results/minimap2/${Name}_merged.bam"
        samtools index "${WD}/results/minimap2/${Name}_merged.bam"

        samtools merge -f "${WD}/results/minimap2/${Name}.bam" \
            "${WD}/results/minimap2/${Name}_PE.bam" \
            "${WD}/results/minimap2/${Name}_merged.bam"
        samtools index "${WD}/results/minimap2/${Name}.bam"

        rm "${WD}/results/minimap2/${Name}_PE.bam"*
        rm "${WD}/results/minimap2/${Name}_merged.bam"*
    fi
done <"${WD}/data/datasets.csv"

# Run mapDamage2 to investigate deamination patterns
conda deactivate
conda activate "${WD}/scripts/mapdamage2"

mkdir -p "${WD}/results/mapDamage" && cd "${WD}/results/mapDamage"

while IFS="," read -r Library Name Age City Country Wolb Type SRA; do
    if [[ "${Library}" != "Library" ]]; then
        echo "Processing library ${Name}"

        mapDamage -i "${WD}/results/minimap2/${Name}.bam" \
            -r "${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz" \
            --rescale \
            --folder="${WD}/results/mapDamage/${Name}"

        convert -density 300 "${WD}/results/mapDamage/${Name}/Stats_out_MCMC_post_pred.pdf" \
            -quality 90 "${WD}/results/mapDamage/${Name}/Stats_out_MCMC_post_pred.png"

        convert -density 300 "${WD}/results/mapDamage/${Name}/Stats_out_MCMC_hist.pdf" \
            -quality 90 "${WD}/results/mapDamage/${Name}/Stats_out_MCMC_hist.png"
    fi
done <"${WD}/data/datasets.csv"

conda deactivate
conda activate "${WD}/scripts/programs"

# QUESTIONS:
# A) How to interpret these plots?
# B) What are the differences between historic and recent samples?

# Investigate mutation frequencies using mpileup and a custom script
gunzip "${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz"
bgzip "${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta"
samtools faidx "${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz"

>"${WD}/results/minimap2/bam.list"
for i in "${WD}/results/minimap2/"*.bam; do
    echo "$i" >>"${WD}/results/minimap2/bam.list"
done

samtools mpileup \
    -B \
    -f "${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz" \
    -r 2L \
    -b "${WD}/results/minimap2/bam.list" >"${WD}/results/minimap2/2L.mpileup"

python "${WD}/scripts/MpileupMutations.py" \
    --mpileup "${WD}/results/minimap2/2L.mpileup" \
    --min-cov 8 \
    --base-quality-threshold 10 \
    --names 18DZ5,19SL15,19SL3,DGRP338,DGRP370,ZI268 >"${WD}/results/minimap2/2L.snps"

less "${WD}/results/minimap2/2L.snps"

# Visualize mutation frequencies in R
mkdir -p "${WD}/results/MutationFreq"

"${WD}/scripts/programs/bin/Rscript" -e "
library(tidyverse)
data <- read.csv('${WD}/results/minimap2/2L.snps', header=FALSE)
colnames(data) <- c('Chrom', 'Pos', 'Library', 'Mutation', 'Freq')
data.condensed <- data %>%
    group_by(Library,Mutation) %>%
    summarise(Count = n()) %>%
    mutate(Freq = Count / sum(Count)) %>%
    separate(Mutation, into = c('Ref', 'Alt'), sep = '>') %>%
    mutate(Ref = toupper(Ref), Alt = toupper(Alt))
ggplot(data.condensed, aes(y = Freq, fill = Library, x=Library)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    facet_grid(Alt~ Ref) +
    labs(title = 'Frequency of Mutations in 2L Region',
         x = 'Library',
         y = 'Frequency') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('${WD}/results/MutationFreq/2L_mutation_frequency.pdf', width = 12, height = 8)
ggsave('${WD}/results/MutationFreq/2L_mutation_frequency.png', width = 12, height = 8)
data <- data %>%
    separate(Mutation, into = c('Ref', 'Alt'), sep = '>') %>%
    mutate(Ref = toupper(Ref), Alt = toupper(Alt)) %>%
    filter(Ref != Alt)
data\$Library <- factor(data\$Library, levels = c('18DZ5', '19SL15', '19SL3', 'DGRP338', 'DGRP370', 'ZI268'))
for (lib in levels(data\$Library)) {
    p <- ggplot(data[data\$Library == lib, ], aes(x = Freq, fill = Library)) +
        geom_histogram(binwidth = 0.05, position = 'dodge') +
        facet_grid(Alt ~ Ref) +
        labs(title = paste('Frequency of Mutations in 2L Region -', lib),
             x = 'Frequency',
             y = 'Count') +
        theme_bw() +
        theme(legend.position = 'none') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0('${WD}/results/MutationFreq/2L_mutation_frequency_histogram_', lib, '.pdf'), p, width = 12, height = 8)
    ggsave(paste0('${WD}/results/MutationFreq/2L_mutation_frequency_histogram_', lib, '.png'), p, width = 12, height = 8)
}
"

# QUESTIONS:
# A) Are these results consistent with mapDamage2?
# B) Are DNA damage signals rare in the datasets?

################################################################################
# 8. Read Depth and Coverage
################################################################################

# Calculate read depths and coverage for each genomic region and library using samtools.

mkdir -p "${WD}/results/ReadDepths" && cd "${WD}/results/ReadDepths"

# Make header line
printf "library\trname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n" \
    >"${WD}/results/ReadDepths/ReadDepths.txt"

# Loop through libraries and calculate summary statistics
for i in "${WD}/results/minimap2/"*.bam; do
    base=$(basename "$i" .bam)
    samtools coverage \
        --reference "${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz" \
        "$i" |
        awk -v base="$base" 'BEGIN{OFS="\t"} NR >1 {print base, $0}' \
            >>"${WD}/results/ReadDepths/ReadDepths.txt"
done

# Plot read depths, coverages, and relative bacterial titer in R
"${WD}/scripts/programs/bin/Rscript" -e "
library(tidyverse)
data <- read.table('${WD}/results/ReadDepths/ReadDepths.txt', header = TRUE, sep = '\t')
data\$library <- factor(data\$library, levels = c('18DZ5', '19SL15', '19SL3', 'DGRP338', 'DGRP370', 'ZI268'))
data\$rname <- factor(data\$rname, levels = c('2L', 'wMel', 'mitochondrion_genome'))
data\$rname <- recode(data\$rname, '2L' = '2L', 'wMel' = 'Wolbachia', 'mitochondrion_genome' = 'Mitochondrion')
ggplot(data, aes(x = rname, y = coverage, fill = library)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    labs(title = 'Coverage Distribution by Region and Library',
         x = 'Region',
         y = 'Coverage') +
    theme_bw() +
    facet_grid(.~library , scales = 'free_y') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
    theme(legend.position = 'bottom')
ggsave('${WD}/results/ReadDepths/coverage_distribution.pdf', width = 8, height = 6)
ggsave('${WD}/results/ReadDepths/coverage_distribution.png', width = 8, height = 6)
ggplot(data, aes(x = rname, y = numreads, fill = library)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    labs(title = 'Read Depth Distribution by Region and Library',
         x = 'Region',
         y = 'Read Depth') +
    theme_bw() +
    facet_grid(.~library , scales = 'free_y') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
    theme(legend.position = 'bottom')
ggsave('${WD}/results/ReadDepths/read_depth_distribution.pdf', width = 8, height = 6)
ggsave('${WD}/results/ReadDepths/read_depth_distribution.png', width = 8, height = 6)
data\$ratio <- ifelse(data\$rname == 'Wolbachia', data\$numreads / data[data\$rname == '2L', ]\$numreads, NA)
ggplot(data, aes(x = library, y = ratio, fill = library)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    labs(title = 'Ratio of Read Depths for Wolbachia and 2L',
         x = 'Library',
         y = 'Ratio of Read Depths') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
    theme(legend.position = 'bottom')
ggsave('${WD}/results/ReadDepths/read_depth_ratio.pdf', width = 8, height = 6)
ggsave('${WD}/results/ReadDepths/read_depth_ratio.png', width = 8, height = 6)
"

# QUESTIONS:
# A) Are there differences in the read depths and coverages between recent and historic samples?
# B) How does the read depth of Wolbachia compare to the nuclear genome (2L)?
# C) What can we infer about the Wolbachia infection status?

################################################################################
# 9. SNP Calling and Phylogenetic Analysis
################################################################################

# Call SNP variants for each genomic region and perform phylogenetic analysis using the rescaled bam files from mapDamage2.

mkdir -p "${WD}/results/SNPs" && cd "${WD}/results/SNPs"

## Use rescaled BAM files for SNP calling
>"${WD}/results/minimap2/scaled_bam.list"
for i in ${WD}/results/mapDamage/*/*.rescaled.bam; do
    echo "$i" >>"${WD}/results/minimap2/scaled_bam.list"
done

# Diploid region (2L)
bcftools mpileup -Ou \
    -f "${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz" \
    -r 2L \
    -a AD,DP \
    -b "${WD}/results/minimap2/scaled_bam.list" |
    bcftools call -mv --ploidy 2 -Ou |
    bcftools view -v snps -Oz -o "${WD}/results/SNPs/2L.diploid.vcf.gz"

python "${WD}/scripts/DiploVCF2Phylip.py" \
    --input "${WD}/results/SNPs/2L.diploid.vcf.gz" \
    --MaxPropGaps 0.1 \
    --MinCov 30 \
    --names 18DZ5,19SL15,19SL3,DGRP338,DGRP370,ZI268 \
    >"${WD}/results/SNPs/2L.phy"

# Haploid region (mitochondrion_genome)
bcftools mpileup -Ou \
    -f "${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz" \
    -r mitochondrion_genome \
    -a AD,DP \
    -b "${WD}/results/minimap2/scaled_bam.list" |
    bcftools call -mv --ploidy 1 -Ou |
    bcftools view -v snps -Oz -o "${WD}/results/SNPs/mito.haploid.vcf.gz"

python "${WD}/scripts/HaploVCF2Phylip.py" \
    --input "${WD}/results/SNPs/mito.haploid.vcf.gz" \
    --MinAlt 1 \
    --MaxPropGaps 0.7 \
    --MinCov 10 \
    --names 18DZ5,19SL15,19SL3,DGRP338,DGRP370,ZI268 \
    >"${WD}/results/SNPs/mito.phy"

# Haploid region (wMel)
bcftools mpileup -Ou \
    -f "${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz" \
    -r wMel \
    -a AD,DP \
    -b "${WD}/results/minimap2/scaled_bam.list" |
    bcftools call -mv --ploidy 1 -Ou |
    bcftools view -v snps -Oz -o "${WD}/results/SNPs/Wolbachia.haploid.vcf.gz"

python "${WD}/scripts/HaploVCF2Phylip.py" \
    --input "${WD}/results/SNPs/Wolbachia.haploid.vcf.gz" \
    --MinAlt 1 \
    --MaxPropGaps 0.5 \
    --MinCov 5 \
    --names 18DZ5,19SL15,19SL3,DGRP338,DGRP370,ZI268 \
    >"${WD}/results/SNPs/Wolbachia.phy"

# Visualize phylogenetic trees in R
"${WD}/scripts/programs/bin/Rscript" -e "
library(ggtree)
library(phangorn)
library(phytools)
library(ape)
library(tidyverse)
library(patchwork)
input_files <- c(
  '${WD}/results/SNPs/2L.phy',
  '${WD}/results/SNPs/mito.phy',
  '${WD}/results/SNPs/Wolbachia.phy'
)
titles <- c('2L 100K Region', 'Mitochondrion', 'Wolbachia')
plots <- list()
for (i in seq_along(input_files)) {
  phylo <- read.phyDat(input_files[i], format = 'phylip')
  tree <- upgma(dist.ml(phylo))
  tree <- midpoint.root(tree)
  Xmax <- max(nodeHeights(tree))
  p <- ggtree(tree, layout = 'roundrect') +
    geom_tiplab(size = 3) +
    ggplot2::xlim(0, Xmax + Xmax/3) +
    ggtitle(titles[i]) +
    theme_tree2() +
    theme_bw()
  plots[[i]] <- p
}
combined_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plot_layout(ncol = 3)
ggsave('${WD}/results/SNPs/combined_phylo_trees.pdf', combined_plot, width = 12, height = 6)
ggsave('${WD}/results/SNPs/combined_phylo_trees.png', combined_plot, width = 12, height = 6)
"

# QUESTIONS:
# A) Are there differences between the nuclear, the mitochondrial and the Wolbachia tree?
# B) How to interpret the long branch in the Wolbachia tree?
# C) What does the tree based on 2L depict?

# Exclude potentially contaminated sample and re-plot phylogenetic trees
python "${WD}/scripts/HaploVCF2Phylip.py" \
    --input "${WD}/results/SNPs/Wolbachia.haploid.vcf.gz" \
    --MinAlt 1 \
    --MaxPropGaps 0.5 \
    --MinCov 5 \
    --exclude 18DZ5,19SL3 \
    --names 18DZ5,19SL15,19SL3,DGRP338,DGRP370,ZI268 \
    >"${WD}/results/SNPs/Wolbachia_no19SL3.phy"

"${WD}/scripts/programs/bin/Rscript" -e "
library(ggtree)
library(phangorn)
library(phytools)
library(ape)
library(tidyverse)
library(patchwork)
input_files <- c(
  '${WD}/results/SNPs/2L.phy',
  '${WD}/results/SNPs/mito.phy',
  '${WD}/results/SNPs/Wolbachia_no19SL3.phy'
)
titles <- c('2L 100K Region', 'Mitochondrion', 'Wolbachia')
plots <- list()
for (i in seq_along(input_files)) {
  phylo <- read.phyDat(input_files[i], format = 'phylip')
  tree <- upgma(dist.ml(phylo))
  Xmax <- max(nodeHeights(tree))
  p <- ggtree(tree, layout = 'roundrect') +
    geom_tiplab(size = 3) +
    ggplot2::xlim(0, Xmax + Xmax/3) +
    ggtitle(titles[i]) +
    theme_tree2() +
    theme_bw()
  plots[[i]] <- p
}
combined_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plot_layout(ncol = 3)
ggsave('${WD}/results/SNPs/combined_phylo_trees_no19SL3.pdf', combined_plot, width = 12, height = 6)
ggsave('${WD}/results/SNPs/combined_phylo_trees_no19SL3.png', combined_plot, width = 12, height = 6)
"

#QUESTIONS:
# A) Why are the mitochondrial and the Wolbachia trees so similar?
# B) What can we infer about the *Wolbachia* variants in the historic sample(s)?

################################################################################
# End of Museomics Workshop 2025 Bioinformatics Pipeline
################################################################################

# Thank you for participating in the Museomics Workshop 2025! I hope you found this bioinformatics pipeline useful for your own research. If you have any questions or feedback, please feel free to reach out to me.
