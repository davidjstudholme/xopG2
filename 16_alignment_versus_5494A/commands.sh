### Download WHRI_5494A reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/035/023/565/GCA_035023565.1_ASM3502356v1/GCA_035023565.1_ASM3502356v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/035/023/565/GCA_035023565.1_ASM3502356v1/GCA_035023565.1_ASM3502356v1_genomic.gff.gz

### Get genomic reads
fastq-dump --split-files --gzip SRR25039548 SRR25039537 SRR25039539 SRR25039538 SRR25039547 SRR25205295

### Ensure that BWA is installed
conda create -n bwa_env
conda activate bwa_env
conda install bioconda::bwa
conda list -n bwa_env > bwa_env_packages.txt
conda env export > bwa_env.yaml

### Align reads against reference genome  
bwa index GCA_035023565.1_ASM3502356v1_genomic.fna.gz
for i in SRR25039548 SRR25039537 SRR25039539 SRR25039538 SRR25039547 SRR25205295
do
    echo bwa mem $i
    bwa mem GCA_035023565.1_ASM3502356v1_genomic.fna.gz $i"_1.fastq.gz" $i"_2.fastq.gz" > $i.versus.WHRI_5494A.sam
    echo convert sam => bam
    samtools view -bS $i.versus.WHRI_5494A.sam -o $i.versus.WHRI_5494A.bam
    echo sort bam
    samtools sort $i.versus.WHRI_5494A.bam -o $i.versus.WHRI_5494A.sorted.bam
    echo rmdup bam
    samtools rmdup $i.versus.WHRI_5494A.sorted.bam $i.versus.WHRI_5494A.sorted.rmdup.bam
    echo index bam
    samtools index $i.versus.WHRI_5494A.sorted.rmdup.bam
    echo Cleaning up intermeidate files
    rm $i.versus.WHRI_5494A.sam $i.versus.WHRI_5494A.bam $i.versus.WHRI_5494A.sorted.bam
    echo done
done



 
