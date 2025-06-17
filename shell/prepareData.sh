WD=/media/inter/mkapun/projects/MuseomicsWorkshop2025

##activate the conda environment
conda activate ${WD}/scripts/PrepareData

## (1) create reference genomes
mkdir -p ${WD}/data/refseq/dmel && cd ${WD}/data/refseq/dmel

# get the latest version of the Drosophila melanogaster genome
wget -O - "https://ftp.flybase.org/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.fa.gz" \
    >${WD}/data/refseq/dmel/dmel-all-chromosome-r6.51.fasta.gz

#subset the genome to 2L:1-100000 and mitochondrion_genome
seqkit subseq --chr 2L -r 1:100000 ${WD}/data/refseq/dmel/dmel-all-chromosome-r6.51.fasta.gz |
    pigz -f >${WD}/data/refseq/dmel/dmel_2L_100K.fasta.gz

seqkit subseq --chr mitochondrion_genome -r 1:100000 ${WD}/data/refseq/dmel/dmel-all-chromosome-r6.51.fasta.gz |
    pigz -f >${WD}/data/refseq/dmel/dmel_mito.fasta.gz

# get the Wolbachia endosymbiont genome variant wMel
wget -qO - "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/584/425/GCF_016584425.1_ASM1658442v1/GCF_016584425.1_ASM1658442v1_genomic.fna.gz" \
    >${WD}/data/refseq/dmel/wMel.fasta.gz

# rename the header of the wMel genome to wMel
pigz -d ${WD}/data/refseq/dmel/wMel.fasta.gz
sed '1s/.*/>wMel/' ${WD}/data/refseq/dmel/wMel.fasta |
    pigz -f >${WD}/data/refseq/dmel/wMel.fasta.gz

# concatenate the 2L, mitochondrion and wMel genomes into one file
cat ${WD}/data/refseq/dmel/dmel_2L_100K.fasta.gz \
    ${WD}/data/refseq/dmel/dmel_mito.fasta.gz \
    ${WD}/data/refseq/dmel/wMel.fasta.gz \
    >${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz

# index the concatenated genome with bwa
bwa index ${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz

## (2) create raw FASTQ files
mkdir -p ${WD}/data/raw_reads && cd ${WD}/data/raw_reads

# (a) loop through the datasets.csv file and subset the reads to 2L:1-100K
while IFS=$"," read -r Library Name Age City Country Wolb Type SRA; do
    # skip the header line
    if [[ "${Library}" != "Library" ]]; then
        echo "Processing library ${Library}"

        # map the reads to the 2L:1-100K genome
        bwa mem \
            -t 150 \
            ${WD}/data/refseq/dmel/dmel_wMel_2L_100K.fasta.gz \
            ${WD}/data/${Type}/${Library}_1.fastq.gz \
            ${WD}/data/${Type}/${Library}_2.fastq.gz |
            samtools view -bS -F 4 - >${WD}/data/${Type}/${Library}.bam

        # isolate the reads that map to 2L:1-100K
        samtools fastq \
            -f 0x2 \
            -F 0x904 \
            -1 ${WD}/data/raw_reads/${Name}_2L_100K_1.fastq.gz \
            -2 ${WD}/data/raw_reads/${Name}_2L_100K_2.fastq.gz \
            ${WD}/data/${Type}/${Library}.bam

    fi

done <${WD}/data/datasets.csv

# trim the reads using fastp
while IFS=$"," read -r Library Name Age City Country Wolb Type SRA; do
    # skip the header line
    if [[ "${Library}" != "Library" ]]; then
        echo "Processing library ${Name}"

        ## trim the reads using fastp
        fastp \
            -i ${WD}/data/raw_reads/${Name}_2L_100K_1.fastq.gz \
            -I ${WD}/data/raw_reads/${Name}_2L_100K_2.fastq.gz \
            -o ${WD}/data/raw_reads/${Name}_2L_100K_1_trimmed.fastq.gz \
            -O ${WD}/data/raw_reads/${Name}_2L_100K_2_trimmed.fastq.gz \
            --merge \
            --merged_out ${WD}/data/raw_reads/${Name}_2L_100K_merged.fastq.gz \
            --length_required 25 \
            --dedup \
            --trim_poly_g \
            --html ${WD}/data/raw_reads/${Name}.html \
            --json ${WD}/data/raw_reads/${Name}.json \
            --detect_adapter_for_pe
    fi

done <${WD}/data/datasets.csv

# (b) subset the reads to 1 million reads
while IFS=$"," read -r Library Name Age City Country Wolb Type SRA; do
    ## skip the header line
    if [[ "${Library}" != "Library" ]]; then
        echo "Processing library ${Library}"

        pigz -dc ${WD}/data/${Type}/${Library}_1.fastq.gz |
            head -n 4000000 |
            pigz -f >${WD}/data/raw_reads/${Name}_1.fastq.gz

        pigz -dc ${WD}/data/${Type}/${Library}_2.fastq.gz |
            head -n 4000000 |
            pigz -f >${WD}/data/raw_reads/${Name}_2.fastq.gz
    fi

done <${WD}/data/datasets.csv

## tar and gzip the raw reads folder
tar -czf ${WD}/data/raw_reads.tar.gz -C ${WD}/data/raw_reads .

## tar and gzip the results folder of mapDamage
tar -czf ${WD}/results/mapDamage.tar.gz -C ${WD}/results/mapDamage .
