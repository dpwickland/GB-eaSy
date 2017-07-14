#!/bin/bash

set -e
. ./GB-eaSy_parameters.txt

##SET UP DIRECTORY STRUCTURE
mkdir -p Intermediate_files/1.Demultiplexed_reads Intermediate_files/2.bam_alignments/ Intermediate_files/3.mpileup/ Intermediate_files/4.Raw_SNPs/ Results/

##DEMULTIPLEX READS
if [ -z $raw_reads_R2 ] #if raw_reads_R2 variable not set in parameters file
	then
		java -jar $GBSX --Demultiplexer -t $num_cores -f1 $raw_reads_R1 -i $barcodes_file -gzip true -ca $adapter_seq -minsl $min_length -o Intermediate_files/1.Demultiplexed_reads; rm Intermediate_files/1.Demultiplexed_reads/*undetermined*
	else
		java -jar $GBSX --Demultiplexer -t $num_cores -f1 $raw_reads_R1 -f2 $raw_reads_R2 -i $barcodes_file -gzip true -ca $adapter_seq -minsl $min_length -o Intermediate_files/1.Demultiplexed_reads; rm Intermediate_files/1.Demultiplexed_reads/*undetermined*
fi

##ALIGN TO REFERENCE
if [ -z $raw_reads_R2 ] #if raw_reads_R2 variable not set in parameters file
	then
		parallel --max-procs $num_cores --keep-order --link "bwa mem  $ref_genome {} | samtools sort -o Intermediate_files/2.bam_alignments/{/.}.sorted_bam; samtools index Intermediate_files/2.bam_alignments/{/.}.sorted_bam" ::: `ls Intermediate_files/1.Demultiplexed_reads/*.R1.fastq.gz` 
	else
		parallel --max-procs $num_cores --keep-order --link "bwa mem  $ref_genome {1} {2} | samtools sort -o Intermediate_files/2.bam_alignments/{1/.}.sorted_bam; samtools index Intermediate_files/2.bam_alignments/{1/.}.sorted_bam" ::: `ls Intermediate_files/1.Demultiplexed_reads/*.R1.fastq.gz` ::: `ls Intermediate_files/1.Demultiplexed_reads/*.R2.fastq.gz`
fi

##CREATE LIST OF SORTED BAM FILES
ls -1 Intermediate_files/2.bam_alignments/*.sorted_bam > Intermediate_files/2.bam_alignments/samples_list.txt

#start GNU parallel to run each region (e.g. chromosome, scaffold) on separate CPU core
parallel  --gnu --max-procs $num_cores --keep-order "\

##GENERATE PILEUP 
bcftools mpileup --regions {} --output-type z --skip-indels --annotate AD,DP --fasta-ref $ref_genome --min-MQ 20 --min-BQ 20  --no-version -b Intermediate_files/2.bam_alignments/samples_list.txt -o Intermediate_files/3.mpileup/mpileup_{}.vcf.gz;\

##CALL SNPS
bcftools call --multiallelic-caller --variants-only --no-version Intermediate_files/3.mpileup/mpileup_{}.vcf.gz | sed -e 's|$(pwd)\/||g' -e 's/Intermediate_files\/2\.bam_alignments\///g' -e  's/\.R.\.fastq.sorted_bam//g'  > Intermediate_files/4.Raw_SNPs/raw_SNPs_{}.vcf;\

" ::: `grep ">" $ref_genome | sed 's/>//g'`

##COMBINE SNPS FROM ALL REGIONS
bcftools concat --no-version `for file in $(ls Intermediate_files/4.Raw_SNPs/*.vcf | sort -V); do echo $file; done;` > Results/all_SNPs_raw.vcf	

##FILTER VCF 
vcftools --vcf Results/all_SNPs_raw.vcf --minDP $min_depth --recode --stdout | awk '/#/ || /[0-9]\/[0-9]/' >Results/all_SNPs_minDP$min_depth.vcf
