# GB-eaSy

## Table of Contents
I. [Introduction](#Introduction)  
II. [Overview](#Overview)  
III. [Software requirements](#Software-requirements)  
IV. [Walkthrough of GB-eaSy pipeline](#Walkthrough-of-GB-eaSy-pipeline)  
V. [Tutorial using sample data](#Tutorial-using-sample-data)  
VI. [Citing GB-eaSy](#Citing-GB-eaSy)  
VII. [License](#License)

## I. Introduction <a name="Introduction"></a>  
GB-eaSy is a GBS bioinformatics pipeline that efficiently incorporates widely used genomics tools, parallelization and automation to increase the accuracy, efficiency and accessibility of GBS analysis. It consists of a Bash shell script that executes several bioinformatics software programs in a Linux environment. This pipeline requires a reference genome and is compatible with both single- and paired-end Illumina reads. It implements the same well-tested tools commonly adopted in whole-genome sequencing.

GB-eaSy is appropriate for users without extensive command-line expertise as well as for experienced bioinformaticians who may choose to modify any step of the script. Before starting the pipeline, the user modifies a parameters file with settings customized for their GBS project (e.g. path to raw sequencer output file, path to barcodes file, number of CPU cores to use). The user then issues a single command to execute the pipeline. 

## II. Overview <a name="Overview"></a>  
### Demultiplex reads and trim adapter sequence
The first step of GB-eaSy uses the software GBSX (Herten et al. 2015) to demultiplex reads and trim adapter sequences based on a user-created barcodes file containing the short barcode sequences that uniquely identify each sample. 
### Align to reference genome
Next, demultiplexed reads are aligned to the reference genome using BWA-MEM (Li 2013). GB-eaSy hastens this alignment step by processing reads files in parallel using GNU parallel (Tange 2001). 
### Call SNPs
After alignment, BCFtools (Li 2011) is used to create a pileup of read bases from which SNPs are then called. This SNP-calling step uses GNU parallel to process each entry in the reference genome file (e.g. each chromosome, each scaffold) on its own CPU core, greatly increasing the efficiency of processing. 
### Filter results 
Finally, the output VCF file is filtered according to a user-specified minimum read depth, using VCFtools (Danecek et al. 2011). 

## III. Software requirements <a name="Software-requirements"></a>  
perl (https://www.perl.org/get.html)  
GNU parallel (version 20170122 or later): http://ftp.gnu.org/gnu/parallel/   
JAVA (version 1.8 or later): https://java.com/en/download/manual.jsp  
GBSX (version 1.3 or later): https://github.com/GenomicsCoreLeuven/GBSX  
BWA (version 0.7.15 or later): https://sourceforge.net/projects/bio-bwa/files/    
zlib (required for bwa): http://www.zlib.net/   
SAMtools (version 1.4 or later): http://www.htslib.org/download/        
BCFtools (version 1.4 or later): http://www.htslib.org/download/  
VCFtools (version 0.1.12b or later): https://sourceforge.net/projects/vcftools/files/     

Versions earlier than those listed above have not been tested and may not work as expected. 
  
## IV. Walkthrough of GB-eaSy pipeline <a name="Walkthrough-of-GB-eaSy-pipeline"></a>  
The following section explains each command that the GB-eaSy.sh script executes. 

### Step 0: Load parameters file
```
#!/bin/bash

set -e
. ./GB-eaSy_parameters.txt
```
**set -e** stops the script if any command produces an error. **GB-eaSy_parameters.txt** contains variables whose values must be customized for a given GBS project. This file must be placed in the same directory as the GB-eaSy.sh script. The following values are specified in the parameters file:  
* Path to reference genome (*REF_GENOME*)
* Path to GBSX.jar file (*GBSX*)
* Path to raw data (*RAW_READS_R1* if single-end reads; *RAW_READS_R1* and *RAW_READS_R2* if paired-end reads)  
* Path to barcodes file (*BARCODES_FILE*)
* Adapter sequence (*ADAPTER_SEQ*)
* Number of CPU cores to use (*NUM_CORES*)
* Minimum read depth to filter SNPs (*MIN_DEPTH*)
* Minimum read length to keep after barcode and adapter trim (*MIN_LENGTH*)  
	
### Step 1: Create directory structure
```
mkdir -p Intermediate_files/1.Demultiplexed_reads Intermediate_files/2.bam_alignments/ Intermediate_files/3.mpileup/ Intermediate_files/4.Raw_SNPs/ Results/
```
In the current directory, two directories are created: **Intermediate_files** (to store demultiplexed fastq files, BAM alignment files, compressed VCF pileup files, and VCF files containing each region's unfiltered SNPs) and **Results** (to store a file containing all identified SNPs and a file containing those SNPs filtered according to read depth).

### Step 2: Demultiplex raw reads
```
if [ -z $RAW_READS_R2 ] #if RAW_READS_R2 variable not set in parameters file, then demultiplex single-end reads
	then
		java -jar $GBSX --Demultiplexer -t $NUM_CORES -f1 $RAW_READS_R1 -i $BARCODES_FILE -gzip true -ca $ADAPTER_SEQ -minsl $MIN_LENGTH -o Intermediate_files/1.Demultiplexed_reads; rm Intermediate_files/1.Demultiplexed_reads/*undetermined*
		
elif [ -n $RAW_READS_R2 ] #if RAW_READS_R2 variable set in parameters file, then demultiplex paired-end reads
	then
		java -jar $GBSX --Demultiplexer -t $NUM_CORES -f1 $RAW_READS_R1 -f2 $RAW_READS_R2 -i $BARCODES_FILE -gzip true -ca $ADAPTER_SEQ -minsl $MIN_LENGTH -o Intermediate_files/1.Demultiplexed_reads; rm Intermediate_files/1.Demultiplexed_reads/*undetermined*
fi
```
**GBSX** is used to separate reads according to the taxa listed in the barcodes file (set by the *BARCODES_FILE* variable in the parameters file). This step has two tracks: If the *RAW_READS_R2* variable is not set in the parameters file, then the first statement is executed to demultiplex single-end reads whose path is set by *RAW_READS_R1*. If the *RAW_READS_R2* (and *RAW_READS_R1*) variable is set, then the second statement is executed to demultiplex paired-end reads. A final step removes the file containing reads with undetermined barcodes (this file would otherwise interfere with later analysis).

The tab-delimited barcodes file used as input for GBSX contains three fields: sample name, barcode and restriction enzyme name. See https://github.com/GenomicsCoreLeuven/GBSX, or the example barcodes file in this repository, for more information.

### Step 3: Align to reference
```
if [ -z $RAW_READS_R2 ] #if RAW_READS_R2 variable not set in parameters file, then align single-end reads to reference genome
	then
		parallel --max-procs $NUM_CORES --keep-order --link "bwa mem  $REF_GENOME {} | samtools sort -o Intermediate_files/2.bam_alignments/{/.}.sorted_bam; samtools index Intermediate_files/2.bam_alignments/{/.}.sorted_bam" ::: `ls Intermediate_files/1.Demultiplexed_reads/*.R1.fastq.gz` 
		
elif [ -n $RAW_READS_R2 ] #if RAW_READS_R2 variable set in parameters file, then align paired-end reads to reference genome
	then
		parallel --max-procs $NUM_CORES --keep-order --link "bwa mem  $REF_GENOME {1} {2} | samtools sort -o Intermediate_files/2.bam_alignments/{1/.}.sorted_bam; samtools index Intermediate_files/2.bam_alignments/{1/.}.sorted_bam" ::: `ls Intermediate_files/1.Demultiplexed_reads/*.R1.fastq.gz` ::: `ls Intermediate_files/1.Demultiplexed_reads/*.R2.fastq.gz`
fi
```
**bwa mem** is used to align reads to the reference genome (whose path is set by the variable *REF_GENOME* in the parameters file), which should be unzipped and then indexed with **bwa index** prior to running GB-eaSy. This step has two tracks: If the *RAW_READS_R2* variable is not set in the parameters file, then the first statement is executed to align single-end reads to the reference genome. If the *RAW_READS_R2* (and *RAW_READS_R1*) variable is set, then the second statement is executed to align paired-end reads to the reference genome. 

After alignment, **samtools sort** outputs sorted BAM alignment files, which are indexed with the command **samtools index**.

To accelerate this process, **GNU parallel** is used to execute the alignment-sort-index command sequence on a separate CPU core for each fastq.gz file contained in the directory Intermediate_files/1.Demultiplexed_reads/.


### Step 4: Create list of sorted BAM files
```
ls -1 Intermediate_files/2.bam_alignments/*.sorted_bam > Intermediate_files/2.bam_alignments/samples_list.txt
```
A list of sorted BAM files is required for the next set of commands.

### Step 5 and 6: Generate pileup and call SNPs	

```
parallel  --gnu --max-procs $NUM_CORES --keep-order "\

bcftools mpileup --regions {} --output-type z --skip-indels --annotate AD,DP --fasta-ref $REF_GENOME --min-MQ 20 --min-BQ 20  --no-version -b Intermediate_files/2.bam_alignments/samples_list.txt -o Intermediate_files/3.mpileup/mpileup_{}.vcf.gz;\

bcftools call --multiallelic-caller --variants-only --no-version Intermediate_files/3.mpileup/mpileup_{}.vcf.gz | sed -e 's|$(pwd)\/||g' -e 's/Intermediate_files\/2\.bam_alignments\///g' -e  's/\.R.\.fastq.sorted_bam//g'  > Intermediate_files/4.Raw_SNPs/raw_SNPs_{}.vcf;\

" ::: `grep ">" $REF_GENOME | cut -d ' ' -f1 | sed 's/>//g'`
```
**bcftools mpileup** is used to generate a pileup of read bases from which **bcftools call** identifies SNPs. Here, **bcftools mpileup** is set to skip indels (and instead look for SNPs), to consider only the bases with a quality score of at least 20 and only the reads with a mapping quality of at least 20, and to output in the compressed vcf format. **bcftools call** is set to use the multiallelic caller algorithm and to output (for a given taxon/sample) only the sites that differ from the reference genome. Two additional steps remove the file extensions from the taxa names in the output VCF files.

To allow parallel processing, **GNU parallel** is set to run the mpileup-call command sequence on a separate core for each region (e.g. each chromosome or each scaffold) listed in the reference genome, as indicated by the command following the three colons at the end of this command block.


### Step 7: Combine SNPS from all regions
```
bcftools concat --no-version `ls -v Intermediate_files/4.Raw_SNPs/*.vcf` > Results/all_SNPs_raw.vcf	
```
The previous step generates one VCF file for each region listed in the reference genome. In step 7, the **bcftools concat** command joins these VCF files together into a single file (all_SNPS_raw.vcf).

### Step 8: Filter VCF
```
vcftools --vcf Results/all_SNPs_raw.vcf --minDP $MIN_DEPTH --recode --stdout | awk '/#/ || /[0-9]\/[0-9]/' >Results/all_SNPs_minDP$MIN_DEPTH.vcf
```

**VCFtools** is used to filter the all_SNPs_raw.vcf file according to the minimum read depth specified in the parameters file. **VCFtools** masks the genotype of any SNP below the minimum read depth with the characters "./." The awk commmand removes any SNPs that do not reach the minimum read depth in any taxon.

## V. Tutorial using sample data <a name="Tutorial-using-sample-data"></a>  
**Step 1: After installing the required software [listed above](#Software-requirements), download this repository** 
```
git clone https://github.com/dpwickland/GB-eaSy.git
```

**Step 2: Set current directory to GB-eaSy directory**
```
cd GB-eaSy/
```

**Step 3: Download the reference genome**  
In this example, we use the soybean reference genome, available at JGI. You may need to set up an account before downloading.

```
#Log in to JGI
#Replace USER_NAME with your user name
#Replace USER_PASSWORD with your password
curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=USER_NAME' --data-urlencode 'password=USER_PASSWORD' -c cookies

#Create directory for reference genome
mkdir reference_genome

#Download reference genome
curl 'http://genome.jgi.doe.gov/ext-api/downloads/get_tape_file?blocking=true&url=/Phytozome/download/_JAMO/56981cf70d87851ee9727dcb/Gmax_275_v2.0.fa.gz' -b cookies > ./reference_genome/Gmax_275_v2.0.fa.gz
```

**Step 4: Unzip and index the reference genome** 

```
#Unzip
gunzip reference_genome/Gmax_275_v2.0.fa.gz

#Index
bwa index reference_genome/Gmax_275_v2.0.fa
```
**Step 5: Modify the GB-easy_parameters.txt file if necessary**  
The variables in the parameter file are set by default to work with this tutorial and the sample data. For this tutorial, the variable *GBSX* points to the GBSX_DW.jar script, which we modified from the original GBSX.jar to include the HindIII cut site, which was not supported initially. Because the sample dataset in this tutorial uses the HindIII enzyme, this modified version of the GBSX script is required here.

**Step 6: Run GB-eaSy**
```
sh ./GB-eaSy.sh
```

Note: For some users' configurations, the [software packages](#Software-requirements) may need to be added to the environmental variable PATH, which tells the bash shell where to find them. One way to do this is to edit the **.bash_profile** configuration file (~/.bash_profile), which exports environmental variables such as the PATH. 
 
For example, if the path to java is /home/user/java/jre1.8.0_131/bin, then the following line should be added to ~/.bash_profile:  
```
export PATH=$PATH:/home/user/java/jre1.8.0_131/bin
```  
If packages are added to the PATH in this manner, then place the following line at the top of the GB-eaSy.sh script:
```
. ~/.bash_profile
```

## VI. Citing GB-eaSy <a name="Citing-GB-eaSy"></a>  
Please cite the following paper in work that uses GB-eaSy:
```
Wickland DP, Battu G, Hudson KA, Diers BW & Hudson ME. (2017). A comparison of genotyping-by-sequencing analysis methods on low-coverage crop datasets shows advantages of a new workflow, GB-eaSy. BMC Bioinformatics, 18:586.
```

## VII. License <a name="License"></a>  
This project is licensed under the terms of the MIT license.


