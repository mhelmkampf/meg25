### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2025                             ###
### 05. Variant calling and SNPs                                             ###
### ======================================================================== ###


### =============================<<< bash >>>===================================


### Establish ssh connection to cluster
ssh user1234@rosa.hpc.uni-oldenburg.de
# When prompted, enter password (invisible)


### Update git repository
cd meg25
git pull


### ALTERNATIVELY: If you receive an error, delete and re-download the repository
# cd
# rm -rf meg25
# git clone https://github.com/mhelmkampf/meg25.git


### Create and navigate to today's working directory
mkdir ~/work/vcal
cd ~/work/vcal



### ============================================================================
### Exercise 1: Assess read quality and trim Illumina reads

module load FastQC
module load cutadapt


### Create link to a pair of Illumina reads as an example
### (H. indigo from Belize, 500 kb each forward and reverse)
ln -s /fs/dss/home/haex1482/share/L18225_indbel_500k_F.fastq L18225_indbel_500k_F.fastq
ln -s /fs/dss/home/haex1482/share/L18225_indbel_500k_R.fastq L18225_indbel_500k_R.fastq


### Run FastQC, a quality control tool for high-throughput sequence reads
### https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
fastqc L18225_indbel_500k_F.fastq


### Open a new terminal window on your local computer and download the HTML results
scp user1234@rosa.hpc.uni-oldenburg.de:/user/user1234/work/vcal/*.html .


### Open HTML results with browser


### Alternatively, print results as text only
# lynx -dump L18225_indbel_500k_F_fastqc.html


### Back on the cluster, remove adapters, low quality bases, and too short reads
### using Cutadapt in paired-end mode
cutadapt -h

cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -q 20,20 \
  -m 120 \
  -o L18225_indbel_trimmed_F.fastq \
  -p L18225_indbel_trimmed_R.fastq \
  L18225_indbel_500k_F.fastq \
  L18225_indbel_500k_R.fastq


### Re-run FastQC on trimmed reads and compare before and after
#>



### ============================================================================
### Exercise 2: Map reads to reference genome

module load BWA
module load SAMtools


### Create link to hamlet reference mitochondrial genome assembly
### (extracted from hybrid assembly of last time)
ln -s /fs/dss/home/haex1482/share/HypPue1_mtg.fas HypPue1_mtg.fas


### Map short reads to reference genome using BWA, a fast and accurate mapping tool
bwa index HypPue1_mtg.fas

bwa mem \
  HypPue1_mtg.fas \
  L18225_indbel_trimmed_F.fastq \
  L18225_indbel_trimmed_R.fastq \
  -t 10 \
  > indbel-mtg_unsorted.sam


### Print first three lines of SAM file
#>


### Convert to BAM, sort, and index
samtools view -Sb indbel-mtg_unsorted.sam | samtools sort -o indbel-mtg.bam
samtools index indbel-mtg.bam


### View alignment
samtools tview indbel-mtg.bam HypPue1_mtg.fas


### Mapping stats
samtools flagstat indbel-mtg.bam


### Optional: calculate average depth
samtools depth -a indbel-mtg.bam | awk '{ sum += $3} END {print "Average depth:", sum/NR }'



### ============================================================================
### Exercise 3: Variant calling with BCFtools

module load BCFtools


### Create links to sorted and indexed BAM files of two additional samples
### (H. gummigutta from Honduras and H. puella from Bocas del Toro, Panama)
ln -s /fs/dss/home/haex1482/share/gumhon-mtg.bam gumhon-mtg.bam
ln -s /fs/dss/home/haex1482/share/pueboc-mtg.bam pueboc-mtg.bam


### Generate genotype likelihoods of all three samples using BCFtools mpileup
### Note on thresholds: -q = mapping quality, -Q = base quality, -d = max read depth
bcftools mpileup \
  -f HypPue1_mtg.fas \
  -q 10 \
  -Q 15 \
  -d 10000 \
  -Ou indbel-mtg.bam gumhon-mtg.bam pueboc-mtg.bam\
  > Joint_mtg.bcf


### Call variants with BCFtools call and generate VCF file
### Note: -m = multiallelic calling algorithm, -v = variant sites only
bcftools call \
  -mv \
  -Ov \
  -o Joint_mtg_snps.vcf \
  Joint_mtg.bcf

bcftools call \
  -m \
  -Ov \
  -o Joint_mtg_all.vcf \
  Joint_mtg.bcf


### To look at the genotypes, find and print last header line (starting with #CHROM)
### plus n following lines (-A n)
grep -A 3 '#CHROM' Joint_mtg_snps.vcf


### Calculate basic summary statistics
bcftools stats Joint_mtg_snps.vcf
# The most interesting part is the 'SN' block (e.g. number of samples, SNPs)
bcftools stats Joint_mtg_snps.vcf | grep 'SN\s'



### ============================================================================
### Solutions:

fastqc L18225_indbel_trimmed_F.fastq

scp user1234@rosa.hpc.uni-oldenburg.de:/user/user1234/work/vcal/*.html .

head -n 3 indbel-mtg_unsorted.sam