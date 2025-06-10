### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2025                             ###
### 08. Selection                                                            ###
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


### Load RStudio module
module load RStudio-Server


### Exectute script to start RStudio
rstudio-start-on-rosa.sh


### Open NEW terminal window and run:
# ssh -N -L 8000: ... (full command see instructions in first terminal window)
# Enter password


### Open http://localhost:8000 in web browser


### Create and navigate to today's working directory
mkdir ~/work/sel
cd ~/work/sel




### ============================================================================
### Exercise 1: 

module load VCFtools


### Create link to today's phased SNP dataset
ln -s /fs/dss/home/haex1482/share/hamlets_LG12_phased.vcf.gz hamlets_LG12_phased.vcf.gz
# Additional parameters: minor allele count 2, thinned to 2 kb, no missing data or indels

# phylo2e_phased.vcf.gz: mac 2


### Check for phasing
zcat hamlets_LG12_phased.vcf.gz | grep -v "##" | head
#> Phased genotypes will use "|" instead of "/"


### Convert phased data to IMPUTE format
vcftools \
  --gzvcf hamlets_LG12_phased.vcf.gz \
  --IMPUTE \
  --out LG12_phased

# This produces:
# - chr1_phased.hap (haplotypes)
# - chr1_phased.legend (like .map)
# - chr1_phased.sample (sample info)


### Rename for selscan
mv LG12_phased.hap chr12.hap
mv LG12_phased.legend chr12.map


### Run selscan (iHS, single population)
selscan --ihs \
        --hap chr1.hap \
        --map chr1.map \
        --out chr1_ihs


### Normalize iHS scores
norm --ihs --bins 20 --files chr1_ihs.ihs.out
# Produces: chr1_ihs.ihs.out.norm


### Print in R


### ===============================<<< R >>>====================================

### Load packages
library(tidyverse)


### Set and check working directory
setwd("~/work/clust")
getwd()


### Read in data for k = 2
admix2 <- read_delim("AncProp_k2_hamlets.tsv", delim = " ", col_names = "Sample")


### Pivot to long format and reformat sample name
long2 <- admix2 %>%
  pivot_longer(cols = starts_with("X"), names_to = "Ancestry", values_to = "Proportion") %>%
  mutate(Name = paste0(str_sub(Sample, -6), "_", str_sub(Sample, 1, 5)))


### Basic plot
p2 <- ggplot(long2, aes(x = Name, y = Proportion, fill = Ancestry)) +   # mandatory: data and variables
    geom_bar(position = "fill", stat = "identity")                      # mandatory: visual representation (geom) 


### Make plot prettier
p2 + scale_fill_manual(values = c("mediumseagreen", "coral")) +         # set colors
  labs(x = NULL, y = NULL, tag = "k = 2") +                             # labels
  theme_minimal() +                                                     # basic style package
  theme(                                                                # specific changes to style (fonts, grid lines etc.)
    text = element_text(color = "grey20"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, color = "gray20", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, color = "gray20"),
    plot.tag = element_text(angle = -90),
    plot.tag.position = c(1.02, 0.6),
    legend.position = "none",
    plot.margin = unit(c(1, 10, 1, 1), "mm")
    )


### Read in and plot ancestry proportions for k = 4 and 6
### Note: you will have to specify additional colors in scale_fill_manual,
### pick one from https://sape.inf.usi.ch/quick-reference/ggplot2/colour
#>



### ============================================================================
### Exercise 2: PCA based on SNP profile

### Install / load packages
install.packages("vcfR")

library(vcfR)
library(adegenet)


### Read VCF file into R
vcf <- read.vcfR("hamlets_LG12_2m2k.vcf")


### Convert from vcfR to genlight object
data <- vcfR2genlight(vcf)


### Principal Component Analysis (PCA)
pca <- glPca(data, nf = 2)   # 2 principal components


### Print principal components to screen
#>


### Convert to tibble, add species information
scores <- as.data.frame(pca$scores) %>%
  rownames_to_column("Sample") %>%
  as_tibble() %>%
  mutate(Species = str_sub(Sample, -6, -4))


### Basic plot
p <- ggplot(data = scores, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 5, alpha = 0.75)


### Make plot prettier
q <- p + scale_color_manual(values = c("mediumseagreen", "orange", "royalblue",
                                   "gray30", "coral", "gray70")) +
  theme_light() +
  theme(
    text = element_text(color = "gray20"),
    panel.grid = element_blank(),
    axis.title.x = element_text(vjust = -1.5)
  )


### Optional: Add labels of interesting samples
labels <- c("54761atlliz", "18267unibel")

(q + geom_text(aes(label = ifelse(Sample %in% labels, Sample, "")),
            size = 3, color = "gray20", vjust = -1)
)


### Optional: Zoom into main cluster
(q + geom_text(aes(label = ifelse(Sample %in% labels, Sample, "")),
            size = 3, color = "gray20", vjust = -1) +
  coord_cartesian(xlim = c(-7, -5), ylim = c(-3, 3))   # zoom into these coordinates
)


### Optional: Plot eigenvalues (proportion of variance explained by each PC)
pca$eig

var <- pca$eig / sum(pca$eig)

barplot(var, main = "Proportion of variance explained", las = 2)



### ============================================================================
### Addendum 1: Plot inbreeding coefficient from session 6

### Read in heterozygosity estimates
het <- read_tsv("../gdiv/Het_hamlets_snp.tsv")


### Add column with species name
hetsp <- het %>%
  mutate(Species = str_sub(INDV, -6, -4))


### Caculate mean chromosome-wide inbreeding coefficient Fis per species
species_fis <- hetsp %>%
  group_by(Species) %>%
  summarise(mean_F = mean(F, na.rm = TRUE)) %>%
  arrange(desc(mean_F))


### Plot Fis
ggplot(species_fis, aes(x = Species, y = mean_F, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("mediumseagreen", "orange", "royalblue", "gray30", "coral", "gray70")) +
  labs(title = NULL, x = NULL, y = "Mean Fis") +
  theme_minimal() +
  theme(
    text = element_text(color = "grey20", size = 14),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )



### ============================================================================
### Addendum 2: Plot nucleotide diversity per species

### Create tibble (data frame) with values (calculated via awk in bash)
species_pi <- tibble(
  Species = c(
    "atlahua",
    "gummigutta",
    "indigo",
    "nigricans",
    "puella",
    "unibel"
  ),
  mean_pi = c(
    0.00799963,
    0.00573394,
    0.00554391,
    0.00550682,
    0.00594079,
    0.00600401
  )
)


### Plot Pi
ggplot(species_pi, aes(x = Species, y = mean_pi, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("mediumseagreen", "orange", "royalblue", "gray30", "coral", "gray70")) +
  labs(title = NULL, x = NULL, y = "Mean Pi") +
  theme_minimal() +
  theme(
    text = element_text(color = "grey20", size = 14),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )



### ============================================================================
### Solutions:

### What is the best k?
cat CV_k2-8.out # CV error (K=2): 0.43992


### Read in data for k = 4 and 6
admix4 <- read_delim("AncProp_k4_hamlets.tsv", delim = " ", col_names = "Sample")
admix6 <- read_delim("AncProp_k6_hamlets.tsv", delim = " ", col_names = "Sample")


### Pivot to long format
long4 <- admix4 %>%
  pivot_longer(cols = starts_with("X"), names_to = "Ancestry", values_to = "Proportion") %>%
    mutate(Name = paste0(str_sub(Sample, -6), "_", str_sub(Sample, 1, 5)))

long6 <- admix6 %>%
  pivot_longer(cols = starts_with("X"), names_to = "Ancestry", values_to = "Proportion") %>%
  mutate(Name = paste0(str_sub(Sample, -6), "_", str_sub(Sample, 1, 5)))


### Basic plot
p4 <- ggplot(long4, aes(x = Name, y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity")

p6 <- ggplot(long6, aes(x = Name, y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity")


### Make plot prettier
p4 + scale_fill_manual(values = c("royalblue", "mediumseagreen", "orange", "coral")) +
  labs(x = NULL, y = NULL, tag = "k = 4") +
  theme_minimal() +
  theme(
    text = element_text(color = "grey20"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, color = "gray20", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, color = "gray20"),
    plot.tag = element_text(angle = -90),
    plot.tag.position = c(1.02, 0.6),
    legend.position = "none",
    plot.margin = unit(c(1, 10, 1, 1), "mm")
    )

p6 + scale_fill_manual(values = c("orange", "gray30", "royalblue",
                                  "violet", "coral", "mediumseagreen", "brown")) +
  labs(x = NULL, y = NULL, tag = "k = 6") +
  theme_minimal() +
  theme(
    text = element_text(color = "grey20"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, color = "gray20", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, color = "gray20"),
    plot.tag = element_text(angle = -90),
    plot.tag.position = c(1.02, 0.6),
    legend.position = "none",
    plot.margin = unit(c(1, 10, 1, 1), "mm")
    )


### Print principal components to screen
pca$scores