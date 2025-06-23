### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2025                             ###
### 08. Detecting selection                                                  ###
### ======================================================================== ###


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bash >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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


### Load RStudio module
module load RStudio-Server


### Exectute script to start RStudio
rstudio-start-on-rosa.sh


### Open NEW terminal window and run:
# ssh -N -L 8000: ... (full command see instructions in first terminal window)
# Enter password


### Open http://localhost:8000 in web browser


### Select "Terminal" to execute bash code from within RStudio


### Create and navigate to today's working directory
mkdir ~/work/sel
cd ~/work/sel



### ============================================================================
### Exercise 1: Scan for Fst outlier loci

module load VCFtools

### Create link to SNP dataset (minor allele count 2, no thinning, without H. altlahua)
ln -s /fs/dss/home/haex1482/share/carib_LG12_snp.vcf.gz carib_LG12_snp.vcf.gz


### Define populations
cp ~/meg25/data/meta/pop* .


# Calculate joint Fst in sliding windows of 50 kb
vcftools \
  --gzvcf carib_LG12_snp.vcf.gz \
  --weir-fst-pop pop_gumhon.txt \
  --weir-fst-pop pop_indbel.txt \
  --weir-fst-pop pop_nigbel.txt \
  --weir-fst-pop pop_puebel.txt \
  --weir-fst-pop pop_unibel.txt \
  --fst-window-step 5000 \
  --fst-window-size 50000 \
  --stdout 1> Fst_lg12_50k.tsv 2> Fst_lg12_50k.log


### Switch to "Console" to execute R code



### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< R >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Load packages
library(tidyverse)


### Set and check working directory
setwd("~/work/sel")
getwd()


### Plot joint Fst, along 50 kb windows (Mantattan plot)
fst_50k <- read_tsv("Fst_lg12_50k.tsv")

(f <- ggplot(data = fst_50k, aes(x = BIN_START, y = WEIGHTED_FST)) +
  geom_point(size = 0.25, alpha = 0.5) +
  labs(x = "Position", y = "Fst") +
  theme_minimal() +
  theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
)


### Find the chromosome-wide, weighted Fst estimate in the VCFtools log file (using bash)
#>


### Assign value to a variable called "chrwide_fst"
#> 


### Add chromosome-wide Fst to plot
f + geom_hline(yintercept = chrwide_fst, color = "blue")


### Identify 99.5% quantile
quantile(fst_50k$WEIGHTED_FST, probs = 0.995)


### Add variable to data table indicating outlier status
out <- fst_50k %>%
  mutate(OUTLIER = case_when(
    WEIGHTED_FST > quantile(WEIGHTED_FST, probs = 0.995) ~ "yes",
    TRUE ~ "no")
  )


### Plot with 99.5% quantile highlighted
(o <- ggplot(data = out, aes(x = BIN_START, y = WEIGHTED_FST, color = OUTLIER)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = chrwide_fst, color = "blue") +
    labs(x = "Window", y = "Fst") +
    scale_color_manual(values = c("gray20", "red")) +
    guides(color = "none") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
)


### Retrieve positions of windows containing Fst peaks
print(out %>% filter(OUTLIER == "yes"), n = nrow(out))


### Genomic region of interest around highest Fst peak
#> 20140001 - 20290000 (150 kb)


### ----------------------------------------------------------------------------
### The following inactive code was used to plot a PCA of the region of interest

### Extract genomic region of interest from VCF file (bash)
# vcftools --gzvcf carib_LG12_snp.vcf.gz \
# --chr LG12 \
# --from-bp 20140001 \
# --to-bp 20290000 \
# --recode \
# --stdout | gzip > fstpeak_region.vcf.gz


### Load packages
# library(vcfR)
# library(adegenet)


### Read VCF file into R
# vcf <- read.vcfR("fstpeak_region.vcf.gz")


### Convert from vcfR to genlight object
# data <- vcfR2genlight(vcf)


### Principal Component Analysis (PCA)
# pca <- glPca(data, nf = 2)


### Convert to tibble, add species information
# scores <- as.data.frame(pca$scores) %>%
#   rownames_to_column("Sample") %>%
#   as_tibble() %>%
#   mutate(Species = str_sub(Sample, -6, -4))


### Plot PCA
# (p <- ggplot(data = scores, aes(x = PC1, y = PC2, color = Species)) +
#   geom_point(size = 5, alpha = 0.75) +
#   scale_color_manual(values = c("orange", "royalblue", "gray30", "coral", "gray70")) +
#   theme_light() +
#   theme(
#     text = element_text(color = "gray20"),
#     panel.grid = element_blank(),
#     axis.title.x = element_text(vjust = -1.5)
#   )
# )


### How do the PCAs of the high Fst region and the whole chromosome differ?
### What trait seems to be affected, pointing to genes related to that trait under selection in the region?



### ============================================================================
### Exercise 2: Identify candidate genes under selection


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< bash >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Create link to the genome annotation file in GFF format
ln -s /fs/dss/home/haex1482/share/HypPue1_annotation.gff HypPue1_annotation.gff


### Extract region of interest from GFF file
awk '$1 == "LG12" && $5 >= 20140001 && $4 <= 20290000' HypPue1_annotation.gff > highfst_region.gff


### What genes are found in this region?


### Find out more by searching for "casz1" at
### NCBI Gene: https://www.ncbi.nlm.nih.gov/gene
### Uniprot: https://www.uniprot.org



### ============================================================================
### Optional: Exercise 3: Extended Haplotype Homozygosity-based test using rehh


### Create link to phased SNP dataset (minor allele count 2, thinned to 2 kb, without H. atlahua)
ln -s /fs/dss/home/haex1482/share/pue_LG12_phased.vcf.gz pue_LG12_phased.vcf.gz
ln -s /fs/dss/home/haex1482/share/uni_LG12_phased.vcf.gz uni_LG12_phased.vcf.gz


### Check for phasing
zcat pue_LG12_phased.vcf.gz | head
#> Phased genotypes use "|" instead of "/"


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< R >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Install / Load packages
install.packages("rehh")

library(rehh)
library(vcfR)


### Read in VCF with phased haplotypes (creates haplohh object)
hap_uni <- data2haplohh(hap_file = "uni_LG12_phased2.vcf.gz",
                        haplotype.in.columns = TRUE,
                        polarize_vcf = TRUE,
                        vcf_reader = "vcfR")


### Scan for EHH (Extended Haplotype Homology) along chromosome
scan_uni <- scan_hh(hap_uni)


### Calculate integrated haplotype score(iHS, standardized across allele frequencies)
ihs_uni <- ihh2ihs(scan_uni)


### Manhattan plot
manhattanplot(ihs_uni)


# Make plot prettier: Remove NAs and filter to iHS results
# ihs_df <- ihs_uni$iHS %>%
#   filter(!is.na(iHS))

# ggplot(ihs_df, aes(x = POSITION, y = iHS)) +
#   geom_point(aes(color = abs(iHS) > 2), size = 1.2) +
#   scale_color_manual(values = c("black", "red"), guide = "none") +
#   geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "red") +
#   labs(title = "iHS along Chromosome (Unicolor)",
#        x = "Genomic Position (bp)",
#        y = "iHS") +
#   theme_minimal()


# Subset the iHS data to the region on LG12
# region <- ihs$ihs %>%
#   filter(POSITION >= 20140001 && POSITION <= 20290000)



### ============================================================================
### Solutions:

### Find the chromosome-wide, weighted Fst estimate in the VCFtools log file (using bash)
#> cat Fst_lg12_50k.log


### Assign value to a variable called "chrwide_fst"
#> chrwide_fst <- 0.059668
