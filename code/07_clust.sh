### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2025                             ###
### 07. Genetic clustering (population structure II)                         ###
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
mkdir ~/work/clust
cd ~/work/clust



### ============================================================================
### Exercise 1: Estimate and plot ancestry proportions (ADMIXTURE)

module load VCFtools
module load PLINK
module load ADMIXTURE


### Create link to today's SNP dataset
ln -s /fs/dss/home/haex1482/share/hamlets_LG12_snp.vcf.gz hamlets_LG12_snp.vcf.gz


### Filter by minor allele count and thin data
vcftools \
    --gzvcf hamlets_LG12_snp.vcf.gz \
    --mac 2 \
    --thin 2000 \
    --recode \
    --stdout | bgzip > hamlets_LG12_2m2k.vcf.gz


### Decompress and remove "LG" from contig names
gzip -cd hamlets_LG12_2m2k.vcf.gz |
    sed 's/LG//g' \
    > hamlets_LG12_2m2k.vcf


### Convert to binary BED format
mkdir bed

plink --vcf hamlets_LG12_2m2k.vcf \
    --make-bed \
    --allow-extra-chr \
    --out bed/hamlets_LG12_2m2k


### Run ADMIXTURE with cross-validation
for k in {2..8}
do
    admixture \
    --cv -j10 \
    bed/hamlets_LG12_2m2k.bed $k > hamlets_LG12_2m2k_k${k}.out
done


### Print CV error to find best k (lowest error)
for k in {2..8}
do
    grep 'CV' hamlets_LG12_2m2k_k${k} \
    >> CV_k2-8.out
done


### Add sample ids to ancestry proportions
for k in {2..8}
do
    paste -d " " ../../meg25/data/meta/hamlet_subset.txt hamlets_LG12_2m2k.${k}.Q |
        sed 's/ $//g' \
        > AncProp_k${k}_hamlets.tsv
done


### To plot the ancestry proportions, we will now switch to R


### Load RStudio module
module load RStudio-Server


### Exectute script to start RStudio
rstudio-start-on-rosa.sh


### Open NEW terminal window and run:
# ssh -N -L 8000: ... (full command see instructions in first terminal window)
# Enter password


### Open http://localhost:8000 in web browser



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
p2 <- ggplot(long2, aes(x = Name, y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity")


### Make plot prettier
p2 + scale_fill_manual(values = c("mediumseagreen", "coral")) +
  labs(x = NULL, y = NULL, tag = "k = 2") +
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


### Read in and plot ancestry proportions for k = 3 and 6
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
vcf


### Convert from vcfR to genlight object
data <- vcfR2genlight(vcf)


### Principal Component Analysis (PCA)
pca <- glPca(data, nf = 2)
pca
pca$scores   # view Principal Components (n = 2)


### Convert to tibble, add species and location information
scores <- as.data.frame(pca$scores) %>%
  rownames_to_column("Sample") %>%
  as_tibble() %>%
  mutate(Species = str_sub(Sample, -6, -4))


### Basic plot
pc <- ggplot(data = scores, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 5, alpha = 0.75)


### Plot PCA
(pc + scale_color_manual(values = c("mediumseagreen", "orange", "royalblue",
                                    "gray30", "coral", "gray70")) +
  theme_light() +
  theme(
    text = element_text(color = "gray20"),
    panel.grid = element_blank(),
    axis.title.x = element_text(vjust = -1.5)
  )
)


### Optional: Plot eigenvalues (proportion of variance explained by each PC)
pca$eig

var <- pca$eig / sum(pca$eig)

barplot(var, main = "Proportion of variance explained", las = 2)


### Optional: zoom into main cluster
(pc + scale_color_manual(values = c("mediumseagreen", "orange", "royalblue", "gray30", "coral", "gray70")) +
  theme_light() +
  theme(
    text = element_text(color = "gray20"),
    panel.grid = element_blank(),
    axis.title.x = element_text(vjust = -1.5)
  ) +
  coord_cartesian(xlim = c(-8, -4), ylim = c(-4, 4))   # zoom in
)



### ============================================================================
### Addon 1: Plot inbreeding coefficient from session 6

### Read in heterozygosity estimates
het <- read_tsv("../gdiv/Het_hamlets_snp.tsv")


### Add column with species name
hetsp <- het %>%
  mutate(Species = str_sub(INDV, -6, -4)
  )

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
### Addon 2: Plot nucleotide diversity per species

### Read in per-site pi estimates

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

### Read in data for k = 3 and 6
admix3 <- read_delim("AncProp_k3_hamlets.tsv", delim = " ", col_names = "Sample")
admix6 <- read_delim("AncProp_k6_hamlets.tsv", delim = " ", col_names = "Sample")


### Pivot to long format
long3 <- admix3 %>%
  pivot_longer(cols = starts_with("X"), names_to = "Ancestry", values_to = "Proportion") %>%
    mutate(Name = paste0(str_sub(Sample, -6), "_", str_sub(Sample, 1, 5)))

long6 <- admix6 %>%
  pivot_longer(cols = starts_with("X"), names_to = "Ancestry", values_to = "Proportion") %>%
  mutate(Name = paste0(str_sub(Sample, -6), "_", str_sub(Sample, 1, 5)))


### Basic plot
p3 <- ggplot(long3, aes(x = Name, y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity")

p6 <- ggplot(long6, aes(x = Name, y = Proportion, fill = Ancestry)) +
    geom_bar(position = "fill", stat = "identity")


### Make plot prettier
p3 + scale_fill_manual(values = c("mediumseagreen", "coral", "royalblue")) +
  labs(x = NULL, y = NULL, tag = "k = 3") +
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