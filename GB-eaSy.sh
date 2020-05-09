#!/bin/bash

PARAMETERS="./GB-eaSy_parameters.txt"

# If a non-default parameters file is supplied, set it as the source of parameters
if [[ -n "$1" ]]; then
	PARAMETERS="$1"
fi

set -e
set -x
. "$PARAMETERS"

##SET UP DIRECTORY STRUCTURE
mkdir -p Intermediate_files/1.Demultiplexed_reads Intermediate_files/2.bam_alignments/ Intermediate_files/3.mpileup/ Intermediate_files/4.Raw_SNPs/ Results/

##DEMULTIPLEX READS
if [ -z $RAW_READS_R2 ] #if RAW_READS_R2 variable not set in parameters file, then demultiplex single-end reads
	then
		java -jar $GBSX --Demultiplexer -t $NUM_CORES -f1 $RAW_READS_R1 -i $BARCODES_FILE -gzip true -ca $ADAPTER_SEQ -minsl $MIN_LENGTH -o Intermediate_files/1.Demultiplexed_reads; rm Intermediate_files/1.Demultiplexed_reads/*undetermined*
		
elif [ -n $RAW_READS_R2 ] #if RAW_READS_R2 variable set in parameters file, then demultiplex paired-end reads
	then
		java -jar $GBSX --Demultiplexer -t $NUM_CORES -f1 $RAW_READS_R1 -f2 $RAW_READS_R2 -i $BARCODES_FILE -gzip true -ca $ADAPTER_SEQ -minsl $MIN_LENGTH -o Intermediate_files/1.Demultiplexed_reads; rm Intermediate_files/1.Demultiplexed_reads/*undetermined*
fi

##ALIGN TO REFERENCE
if [ -z $RAW_READS_R2 ] #if RAW_READS_R2 variable not set in parameters file, then align single-end reads to reference genome
	then
		parallel --max-procs $NUM_CORES --keep-order --link "bwa mem  $REF_GENOME {} | samtools sort -o Intermediate_files/2.bam_alignments/{/.}.sorted_bam; samtools index Intermediate_files/2.bam_alignments/{/.}.sorted_bam" ::: `ls Intermediate_files/1.Demultiplexed_reads/*.R1.fastq.gz` 
		
elif [ -n $RAW_READS_R2 ] #if RAW_READS_R2 variable set in parameters file, then align paired-end reads to reference genome
	then
		parallel --max-procs $NUM_CORES --keep-order --link "bwa mem  $REF_GENOME {1} {2} | samtools sort -o Intermediate_files/2.bam_alignments/{1/.}.sorted_bam; samtools index Intermediate_files/2.bam_alignments/{1/.}.sorted_bam" ::: `ls Intermediate_files/1.Demultiplexed_reads/*.R1.fastq.gz` ::: `ls Intermediate_files/1.Demultiplexed_reads/*.R2.fastq.gz`
fi

##CREATE LIST OF SORTED BAM FILES
ls -1 Intermediate_files/2.bam_alignments/*.sorted_bam > Intermediate_files/2.bam_alignments/samples_list.txt

#start GNU parallel to run each region (e.g. chromosome, scaffold) on separate CPU core
parallel  --gnu --max-procs $NUM_CORES --keep-order "\

##GENERATE PILEUP 
bcftools mpileup --regions {} --output-type z --skip-indels --annotate AD,DP --fasta-ref $REF_GENOME --min-MQ 20 --min-BQ 20  --no-version -b Intermediate_files/2.bam_alignments/samples_list.txt -o Intermediate_files/3.mpileup/mpileup_{}.vcf.gz;\

##CALL SNPS
bcftools call --multiallelic-caller --variants-only --no-version Intermediate_files/3.mpileup/mpileup_{}.vcf.gz | sed -e 's|$(pwd)\/||g' -e 's/Intermediate_files\/2\.bam_alignments\///g' -e  's/\.R.\.fastq.sorted_bam//g'  > Intermediate_files/4.Raw_SNPs/raw_SNPs_{}.vcf;\

" ::: `grep ">" $REF_GENOME | cut -d ' ' -f1 | sed 's/>//g'`

##COMBINE SNPS FROM ALL REGIONS
bcftools concat --no-version `ls -v Intermediate_files/4.Raw_SNPs/*.vcf` > Results/all_SNPs_raw.vcf	

##FILTER VCF 
vcftools --vcf Results/all_SNPs_raw.vcf --minDP $MIN_DEPTH --recode --stdout | awk '/#/ || /[0-9]\/[0-9]/' >Results/all_SNPs_minDP$MIN_DEPTH.vcf
