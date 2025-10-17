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

# Set Working directory
WD="<your/path/to/Workshop_VI_Museomics>"
# WD="/home/mkapun/github/Workshop_VI_Museomics"

# Go to working directory
cd ${WD}

# copy scripts and programs 
cp -r /media/scratch/Museomics_WS_stuff/scripts .

################################################################################
# 2. Copy and Decompress Raw Input FASTQ Data
################################################################################

## Copy and extract raw sequencing reads.
mkdir -p ${WD}/{data,results,QSUB}
mkdir -p ${WD}/data/raw_reads

# Untar the compressed folder within ${WD}/data
cp /media/scratch/Museomics_WS_stuff/data/raw_reads.tar.gz data/
tar -xzf ${WD}/data/raw_reads.tar.gz -C ${WD}/data/raw_reads

## Optionally download mapDamage2 results
mkdir -p ${WD}/results/mapDamage
cp /media/scratch/Museomics_WS_stuff/data/mapDamage2.tar.gz data/
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
echo """Library,Name,Age,City,Country,Wolb,Type,SRA
DGRP370,DGRP370,2003,Raleigh,USA,wMel,recent,SRR834539
DGRP338,DGRP338,2003,Raleigh,USA,wMelCS,recent,SRR834513
ZI268,ZI268,2003,Ziawonga,Zambia,wMel,recent,SRR189425
HG_15,19SL15,1933,Lund,Sweden,unknown,historic,SRR23876580
HG0027,19SL3,1933,Lund,Sweden,unknown,historic,SRR23876574
HG0029,18DZ5,1899,Zealand,Denmark,unknown,historic,SRR23876565""" >${WD}/data/datasets.csv


################################################################################
# 4. Preview Raw FASTQ Data
################################################################################

# Preview the first 10 lines of the forward read file for two datasets.
less ${WD}/data/raw_reads/ZI268_1.qfastq.gz
less ${WD}/data/raw_reads/18DZ5_1.fastq.gz

# QUESTIONS:
# A) What is the structure of the datasets?
# B) How do the two files differ?

################################################################################
# 5. Trim Reads with fastp
################################################################################

# Trim and filter reads using fastp. Focus on a subset of 1,000,000 reads from each dataset.

######## Load dependencies #######

# First load fastp into your environment
mkdir -p ${WD}/data/trimmed_reads && cd ${WD}/data/trimmed_reads

# Loop through each library in datasets.csv and trim reads
while IFS=$"," read -r Library Name Age City Country Wolb Type SRA; do
    if [[ ${Library} != "Library" ]]; then

        echo """
        #!/bin/sh

        ## name of Job
        #PBS -N FASTP_${Name}

        ## Redirect output stream to this file.
        #PBS -o ${WD}/data/trimmed_reads/${Name}_fastp.log

        ## Stream Standard Output AND Standard Error to outputfile (see above)
        #PBS -j oe

        ## Select a maximum of 5 cores and 50gb of RAM
        #PBS -l select=1:ncpus=5:mem=20gb

        ##### load dependencies 
        source /opt/anaconda3/etc/profile.d/conda.sh
        conda activate ${WD}/scripts/programs

        fastp \
            -i ${WD}/data/raw_reads/${Name}_1.fastq.gz \
            -I ${WD}/data/raw_reads/${Name}_2.fastq.gz \
            -o ${WD}/data/trimmed_reads/${Name}_1_trimmed.fastq.gz \
            -O ${WD}/data/trimmed_reads/${Name}_2_trimmed.fastq.gz \
            --merge \
            --merged_out ${WD}/data/trimmed_reads/${Name}_merged.fastq.gz \
            --length_required 25 \
            --thread 5 \
            --dedup \
            --trim_poly_g \
            --html ${WD}/data/trimmed_reads/${Name}.html \
            --json ${WD}/data/trimmed_reads/${Name}.json \
            --detect_adapter_for_pe
        """ >${WD}/QSUB/fastp_${Name}.sh
        
        # send to OpenPBS
        qsub ${WD}/QSUB/fastp_${Name}.sh
    fi
    
done <${WD}/data/datasets.csv

# View HTML output files
firefox ${WD}/data/trimmed_reads/18DZ5.html
firefox ${WD}/data/trimmed_reads/ZI268.html

# QUESTIONS:
# A) What do the parameters mean?
# B) How do the read length distributions differ between the datasets?

################################################################################
# 6. Testing for Eukaryotic Contamination
################################################################################

# Download mitochondrial genomes for contamination screening and build a BLAST database.
mkdir -p ${WD}/data/refseq/contamination && cd ${WD}/data/refseq/contamination

# Download and rename mitochondrial genomes
wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=PP230540.1&db=nuccore&report=fasta&format=text" -O Gryllus_bimaculatus_mito.fasta
sed -i '1s/.*/>Gryllus_bimaculatus/' Gryllus_bimaculatus_mito.fasta

wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_012920.1&db=nuccore&report=fasta&format=text" -O Homo_sapiens_mito.fasta
sed -i '1s/.*/>Homo_sapiens/' Homo_sapiens_mito.fasta

wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_024511.2&db=nuccore&report=fasta&format=text" -O Drosophila_melanogaster_mito.fasta
sed -i '1s/.*/>Drosophila_melanogaster/' Drosophila_melanogaster_mito.fasta

wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=PV608434.1&db=nuccore&report=fasta&format=text" -O Anthrenus_verbasci_mito.fasta
sed -i '1s/.*/>Anthrenus_verbasci/' Anthrenus_verbasci_mito.fasta

# Concatenate all files into one
cat *_mito.fasta >mitochondrial.fasta

# make BLAST database
echo """
#!/bin/sh

## name of Job
#PBS -N MakeDB

## Redirect output stream to this file.
#PBS -o ${WD}/data/refseq/contamination/MakeDB.log

## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe

## Select a maximum of 5 cores and 50gb of RAM
#PBS -l select=1:ncpus=5:mem=20gb

##### load dependencies 
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate ${WD}/scripts/programs

# Build nucleotide BLAST database
makeblastdb \
    -in ${WD}/data/refseq/contamination/mitochondrial.fasta \
    -dbtype nucl \
    -out ${WD}/data/refseq/contamination/mitoDB

""" >${WD}/QSUB/makeblastdb.sh

#send to OpenPBS
qsub ${WD}/QSUB/makeblastdb.sh

# Convert trimmed FASTQ files to FASTA and run BLAST against the mitochondrial database.

# Prepare results folder and output file
mkdir -p ${WD}/results/contamination
>${WD}/results/contamination/BLAST.tsv

# Loop through all files in dataset, convert FASTQ to FASTA, then BLAST
while IFS="," read -r Library Name Age City Country Wolb Type SRA; do
    if [[ ${Library} != "Library" ]]; then
        gunzip -c ${WD}/data/trimmed_reads/${Name}_1_trimmed.fastq.gz |
            awk 'NR % 4 == 1 {print ">" substr($0, 2)} NR % 4 == 2 {print $0}' >${WD}/data/trimmed_reads/${Name}_trimmed.fasta

        gunzip -c ${WD}/data/trimmed_reads/${Name}_merged.fastq.gz |
            awk 'NR % 4 == 1 {print ">" substr($0, 2)} NR % 4 == 2 {print $0}' >>${WD}/data/trimmed_reads/${Name}_trimmed.fasta

        echo """
        #!/bin/sh

        ## name of Job
        #PBS -N BLAST_${Name}

        ## Redirect output stream to this file.
        #PBS -o ${WD}/data/refseq/contamination/${NAME}_blast.log

        ## Stream Standard Output AND Standard Error to outputfile (see above)
        #PBS -j oe

        ## Select a maximum of 5 cores and 50gb of RAM
        #PBS -l select=1:ncpus=5:mem=20gb

        ##### load dependencies 
        source /opt/anaconda3/etc/profile.d/conda.sh
        conda activate ${WD}/scripts/programs

        blastn \
            -num_threads 4 \
            -evalue 1e-10 \
            -max_target_seqs 1 \
            -outfmt '6 qseqid sseqid slen qlen pident length mismatch evalue bitscore' \
            -db ${WD}/data/refseq/contamination/mitoDB \
            -query ${WD}/data/trimmed_reads/${Name}_trimmed.fasta |
            awk -v Name='${Name}' '\$5>90 {print Name\"\t\"\$0}' >>${WD}/results/contamination/BLAST.tsv
        
        """ >${WD}/QSUB/blast_${Name}.sh

        # send to OpenPBS
        qsub -W block=true ${WD}/QSUB/blast_${Name}.sh
    fi
done <${WD}/data/datasets.csv

# View tabular output file
less ${WD}/results/contamination/BLAST.tsv

# QUESTIONS:
# A) What do the columns mean?
# B) Do you see non-Drosophila hits?

# Plot the proportion of reads mapped to mitochondrial genomes in R
${WD}/scripts/programs/bin/Rscript -e "
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
${WD}/scripts/programs/bin/Rscript -e "
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
# 7. ECMSD Analysis
################################################################################

echo """
#!/bin/sh

## name of Job
#PBS -N ECMSD_19SL3

## Redirect output stream to this file.
#PBS -o ${WD}/data/refseq/contamination/19SL3_ecmsd.log

## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe

## Select a maximum of 5 cores and 50gb of RAM
#PBS -l select=1:ncpus=5:mem=20gb

bash /media/inter/pipelines/ECMSD/shell/ECMSD.sh \
    --fwd ${WD}/data/trimmed_reads/19SL3_merged.fastq.gz \
    --out ${WD}/results/contamination/ECMSD/19SL3 \
    --threads 5 \
    --Binsize 1000 \
    --RMUS-threshold 0.15 \
    --mapping_quality 20 \
    --taxonomic-hierarchy genus \
    --force

""" >${WD}/QSUB/ecmsd_19SL3.sh

#send to OpenPBS
qsub ${WD}/QSUB/ecmsd_19SL3.sh

################################################################################
# 8. Testing for DNA Degradation
################################################################################

# Map reads to reference genomes using minimap2 and analyze DNA damage patterns with mapDamage2.
mkdir -p ${WD}/results/minimap2 && cd ${WD}/results/minimap2
while IFS="," read -r Library Name Age City Country Wolb Type SRA; do
    if [[ ${Library} != "Library" ]]; then
        
        echo """
        #!/bin/sh

        ## name of Job
        #PBS -N Map_${Name}

        ## Redirect output stream to this file.
        #PBS -o ${WD}/results/minimap2/${Name}_map.log

        ## Stream Standard Output AND Standard Error to outputfile (see above)
        #PBS -j oe

        ## Select a maximum of 5 cores and 50gb of RAM
        #PBS -l select=1:ncpus=5:mem=20gb

        ##### load dependencies
        source /opt/anaconda3/etc/profile.d/conda.sh
        conda activate ${WD}/scripts/programs

        ## map reads to reference genome
        minimap2 -ax sr --secondary=no -t 5 \
            ${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz \
            ${WD}/data/trimmed_reads/${Name}_1_trimmed.fastq.gz \
            ${WD}/data/trimmed_reads/${Name}_2_trimmed.fastq.gz |
            samtools view -bS -F 4 - |
            samtools sort -o ${WD}/results/minimap2/${Name}_PE.bam
        samtools index ${WD}/results/minimap2/${Name}_PE.bam

        minimap2 -ax sr --secondary=no -t 5 \
            ${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz \
            ${WD}/data/trimmed_reads/${Name}_merged.fastq.gz |
            samtools view -bS -F 4 - |
            samtools sort -o ${WD}/results/minimap2/${Name}_merged.bam
        samtools index ${WD}/results/minimap2/${Name}_merged.bam

        samtools merge -f ${WD}/results/minimap2/${Name}.bam \
            ${WD}/results/minimap2/${Name}_PE.bam \
            ${WD}/results/minimap2/${Name}_merged.bam
        samtools index ${WD}/results/minimap2/${Name}.bam

        rm ${WD}/results/minimap2/${Name}_PE.bam*
        rm ${WD}/results/minimap2/${Name}_merged.bam*

        """ >${WD}/QSUB/map_${Name}.sh

        # send to OpenPBS
        qsub ${WD}/QSUB/map_${Name}.sh

    fi
done <${WD}/data/datasets.csv

# Run mapDamage2 to investigate deamination patterns
mkdir -p ${WD}/results/mapDamage && cd ${WD}/results/mapDamage
while IFS="," read -r Library Name Age City Country Wolb Type SRA; do
    if [[ ${Library} != "Library" ]]; then

        echo """
        #!/bin/sh

        ## name of Job
        #PBS -N MapDamage_${Name}

        ## Redirect output stream to this file.
        #PBS -o ${WD}/results/mapDamage/${Name}_mapDamage.log

        ## Stream Standard Output AND Standard Error to outputfile (see above)
        #PBS -j oe

        ## Select a maximum of 5 cores and 50gb of RAM
        #PBS -l select=1:ncpus=5:mem=20gb

        ##### load dependencies
        source /opt/anaconda3/etc/profile.d/conda.sh
        conda activate ${WD}/scripts/mapdamage2

        mapDamage -i ${WD}/results/minimap2/${Name}.bam \
            -r ${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz \
            --rescale \
            --folder=${WD}/results/mapDamage/${Name}

        ## convert PDFs to PNGs (only works on Linux systems)
        for pdf in ${WD}/results/mapDamage/${Name}/*.pdf; do
            png=\${pdf%.pdf}.png
            convert -density 300 \$pdf -quality 90 \$png
        done

        """ >${WD}/QSUB/mapDamage_${Name}.sh

        # send to OpenPBS
        qsub ${WD}/QSUB/mapDamage_${Name}.sh
    fi
done <${WD}/data/datasets.csv

# Let's summarize the results in a table across all libraries based on the Stats_out_MCMC_iter_summ_stat.csv files

# make output folder
mkdir -p ${WD}/results/mapDamage && cd ${WD}/results/mapDamage
# make empty output file
echo "Name Theta DeltaD DeltaS Lambda" >${WD}/results/mapDamage/MapDamage_summary.txt

# loop through libraries in datasets.csv
while IFS="," read -r Library Name Age City Country Wolb Type SRA; do
    if [[ ${Library} != "Library" ]]; then
        echo "Processing library ${Name}"
        awk -F "," -v V=${Name} '/Mean/ {print V, $2, $3, $4, $5}' ${WD}/results/mapDamage/${Name}/Stats_out_MCMC_iter_summ_stat.csv \
            >>${WD}/results/mapDamage/MapDamage_summary.txt
    fi
done <${WD}/data/datasets.csv

# Now plot the mean values across all libraries
${WD}/scripts/programs/bin/Rscript -e "
library(tidyverse)

## Load the data
data <- read.table('${WD}/results/mapDamage/MapDamage_summary.txt', header = TRUE, sep = ' ')
data.long<-reshape(data,
    direction='long',
    varying=list(colnames(data)[2:ncol(data)]),
    timevar='Parameter',
    times=colnames(data)[2:ncol(data)])

# plot the coverage distribution by region and library
PLOT<-ggplot(data.long, aes(x = Name, y = Theta, fill = Name)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    labs(title = 'Briggs Paramters',
         x = 'Name',
         y = 'Value') +
    theme_bw() +
    facet_wrap(.~Parameter , scales = 'free_y') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot as PDF and PNG
ggsave('${WD}/results/mapDamage/MapDamage_summary.pdf',PLOT, width = 8, height = 4)
ggsave('${WD}/results/mapDamage/MapDamage_summary.png',PLOT, width = 8, height = 4)
"

# QUESTIONS:
# A) How to interpret these plots?
# B) What are the differences between historic and recent samples?

################################################################################
# 9. Read Depth and Coverage
################################################################################

# Calculate read depths and coverage for each genomic region and library using samtools.
mkdir -p ${WD}/results/ReadDepths && cd ${WD}/results/ReadDepths

# Make header line
printf "library\trname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n" \
    >${WD}/results/ReadDepths/ReadDepths.txt

# Loop through libraries and calculate summary statistics
for i in ${WD}/results/minimap2/*.bam; do
    base=$(basename $i .bam)

    echo """
    #!/bin/sh

    ## name of Job
    #PBS -N coverage_${base}

    ## Redirect output stream to this file.
    #PBS -o ${WD}/results/ReadDepths/${base}_coverage.log

    ## Stream Standard Output AND Standard Error to outputfile (see above)
    #PBS -j oe

    ## Select a maximum of 5 cores and 50gb of RAM
    #PBS -l select=1:ncpus=5:mem=20gb

    ##### load dependencies
    source /opt/anaconda3/etc/profile.d/conda.sh
    conda activate ${WD}/scripts/programs

    samtools coverage \
        --reference ${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz \
        $i |
        awk -v base=\"$base\" 'BEGIN{OFS=\"\t\"} NR >1 {print base, \$0}' \
            >>${WD}/results/ReadDepths/ReadDepths.txt
    """ >${WD}/QSUB/coverage_${base}.sh

    # send to OpenPBS
    qsub -W block=true ${WD}/QSUB/coverage_${base}.sh

done

# Plot read depths, coverages, and relative bacterial titer in R
"${WD}/scripts/programs/bin/Rscript" -e "
library(tidyverse)
data <- read.table('${WD}/results/ReadDepths/ReadDepths.txt', header = TRUE, sep = '\t')
data\$library <- factor(data\$library, levels = c('18DZ5', '19SL15', '19SL3', 'DGRP338', 'DGRP370', 'ZI268'))
data\$rname <- factor(data\$rname, levels = c('2L', 'wMel', 'mitochondrion_genome'))
data\$rname <- recode(data\$rname, '2L' = '2L', 'wMel' = 'Wolbachia', 'mitochondrion_genome' = 'Mitochondrion')
PLOT<-ggplot(data, aes(x = rname, y = coverage, fill = library)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    labs(title = 'Coverage Distribution by Region and Library',
         x = 'Region',
         y = 'Coverage') +
    theme_bw() +
    facet_grid(.~library , scales = 'free_y') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
    theme(legend.position = 'bottom')
ggsave('${WD}/results/ReadDepths/coverage_distribution.pdf',PLOT, width = 8, height = 6)
ggsave('${WD}/results/ReadDepths/coverage_distribution.png',PLOT, width = 8, height = 6)
PLOT<-ggplot(data, aes(x = rname, y = numreads, fill = library)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    labs(title = 'Read Depth Distribution by Region and Library',
         x = 'Region',
         y = 'Read Depth') +
    theme_bw() +
    facet_grid(.~library , scales = 'free_y') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
    theme(legend.position = 'bottom')
ggsave('${WD}/results/ReadDepths/read_depth_distribution.pdf', PLOT,width = 8, height = 6)
ggsave('${WD}/results/ReadDepths/read_depth_distribution.png', PLOT,width = 8, height = 6)
data\$ratio <- ifelse(data\$rname == 'Wolbachia', data\$numreads / data[data\$rname == '2L', ]\$numreads, NA)
PLOT<-ggplot(data, aes(x = library, y = ratio, fill = library)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    labs(title = 'Ratio of Read Depths for Wolbachia and 2L',
         x = 'Library',
         y = 'Ratio of Read Depths') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
    theme(legend.position = 'bottom')
ggsave('${WD}/results/ReadDepths/read_depth_ratio.pdf',PLOT, width = 8, height = 6)
ggsave('${WD}/results/ReadDepths/read_depth_ratio.png',PLOT, width = 8, height = 6)
"

# QUESTIONS:
# A) Are there differences in the read depths and coverages between recent and historic samples?
# B) How does the read depth of Wolbachia compare to the nuclear genome (2L)?
# C) What can we infer about the Wolbachia infection status?

################################################################################
# 9. SNP Calling and Phylogenetic Analysis
################################################################################

# Call SNP variants for each genomic region and perform phylogenetic analysis using the rescaled bam files from mapDamage2.
mkdir -p ${WD}/results/SNPs && cd ${WD}/results/SNPs

##### load dependencies
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate ${WD}/scripts/programs

## Use rescaled BAM files for SNP calling
>${WD}/results/minimap2/scaled_bam.list
for i in ${WD}/results/mapDamage/*/*.rescaled.bam; do
    echo "$i" >>${WD}/results/minimap2/scaled_bam.list

    ## make index for each rescaled BAM file
    samtools index "$i"
done

conda deactivate

## do the actual calling
echo """
#!/bin/sh

## name of Job
#PBS -N SNPcalling

## Redirect output stream to this file.
#PBS -o ${WD}/results/SNPs/SNPcalling.log

## Stream Standard Output AND Standard Error to outputfile (see above)
#PBS -j oe

## Select a maximum of 5 cores and 50gb of RAM
#PBS -l select=1:ncpus=5:mem=20gb

##### load dependencies
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate ${WD}/scripts/programs

# Diploid region (2L)
bcftools mpileup -Ou \
    -f ${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz \
    -r 2L \
    -a AD,DP \
    -b ${WD}/results/minimap2/scaled_bam.list |
    bcftools call -mv --ploidy 2 -Ou |
    bcftools view -v snps -Oz -o ${WD}/results/SNPs/2L.diploid.vcf.gz

python3 ${WD}/scripts/DiploVCF2Phylip.py \
    --input ${WD}/results/SNPs/2L.diploid.vcf.gz \
    --MaxPropGaps 0.1 \
    --MinCov 30 \
    >${WD}/results/SNPs/2L.phy

# Haploid region (mitochondrion_genome)
bcftools mpileup -Ou \
    -f ${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz \
    -r mitochondrion_genome \
    -a AD,DP \
    -b ${WD}/results/minimap2/scaled_bam.list |
    bcftools call -mv --ploidy 1 -Ou |
    bcftools view -v snps -Oz -o ${WD}/results/SNPs/mito.haploid.vcf.gz

python3 ${WD}/scripts/HaploVCF2Phylip.py \
    --input ${WD}/results/SNPs/mito.haploid.vcf.gz \
    --MinAlt 1 \
    --MaxPropGaps 0.7 \
    --MinCov 10 \
    >${WD}/results/SNPs/mito.phy

# Haploid region (wMel)
bcftools mpileup -Ou \
    -f ${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz \
    -r wMel \
    -a AD,DP \
    -b ${WD}/results/minimap2/scaled_bam.list |
    bcftools call -mv --ploidy 1 -Ou |
    bcftools view -v snps -Oz -o ${WD}/results/SNPs/Wolbachia.haploid.vcf.gz

python3 ${WD}/scripts/HaploVCF2Phylip.py \
    --input ${WD}/results/SNPs/Wolbachia.haploid.vcf.gz \
    --MinAlt 1 \
    --MaxPropGaps 0.5 \
    --MinCov 5 \
    >${WD}/results/SNPs/Wolbachia.phy

""" >${WD}/QSUB/snp_calling.sh

# send to OpenPBS
qsub ${WD}/QSUB/snp_calling.sh

# Visualize phylogenetic trees in R
${WD}/scripts/programs/bin/Rscript -e "
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
python3 ${WD}/scripts/HaploVCF2Phylip.py \
    --input ${WD}/results/SNPs/Wolbachia.haploid.vcf.gz \
    --MinAlt 1 \
    --MaxPropGaps 0.5 \
    --MinCov 5 \
    --exclude 18DZ5.rescaled,19SL3.rescaled \
    >${WD}/results/SNPs/Wolbachia_no19SL3.phy

${WD}/scripts/programs/bin/Rscript -e "
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

# QUESTIONS:
# A) Why are the mitochondrial and the Wolbachia trees so similar?
# B) What can we infer about the *Wolbachia* variants in the historic sample(s)?

################################################################################
# End of Museomics Workshop 2025 Bioinformatics Pipeline
################################################################################

# Thank you for participating in the Museomics Workshop 2025! I hope you found this bioinformatics pipeline useful for your own research. If you have any questions or feedback, please feel free to reach out to me.
