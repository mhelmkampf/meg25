### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2025                             ###
### 06. Population genomics and genetic diversity                            ###
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
mkdir ~/work/gdiv
cd ~/work/gdiv



### ============================================================================
### Exercise 1: Summarize and filter VCF files

module load BCFtools
module load VCFtools
# Manual: https://samtools.github.io/bcftools/bcftools.html


### Create link to today's SNP dataset
ln -s /fs/dss/home/haex1482/share/hamlets_LG12_snp.vcf.gz hamlets_LG12_snp.vcf.gz


### How many sites (SNPs) does the example dataset have?


### VCFtools is a standard software to filter and run calculations on genetic
### variation data in VCF format: https://vcftools.github.io/man_latest.html


### Use the following code block to filter and write data to a new file
### (insert filtering parameters in lieu of ...)

vcftools \
    --gzvcf ../../data/snps_hamlets_lg12.vcf.gz \
    ... \
    --recode \
    --stdout | bgzip > snps_hamlets_filtered.vcf.gz


### Referring to the manual linked above, include only sites with a "Minor Allele 
### Count" (mac) of at least 2


### Re-calculate the number of sites to confirm filtering was sucessful, and
### spot check visually that there are indeed at least 2 counts of each allele


### ============================================================================
### Exercise 2: Filter by linkage disequilibrium

### Referring to the linkage decay plot, how many basepairs (bp) apart should the 
### sites be from each other to avoid linkage disequilibrium? Use vcftools --thin


### Calculate the correlation coefficient of linkage disequilibrium in the first
### 50000 bp of both datasets (before and after thinning)
vcftools --gzvcf snps_hamlets_filtered.vcf.gz \
  --chr LG12 \
  --from-bp 1 \
  --to-bp 100000 \
  --hap-r2 \
  --stdout > R2_before_thinning.tsv

vcftools --gzvcf snps_hamlets_2kb.vcf.gz \
  --chr LG12 \
  --from-bp 1 \
  --to-bp 100000 \
  --hap-r2 \
  --stdout > R2_after_thinning.tsv


### Load RStudio module
module load RStudio-Server


### Exectute script to start RStudio
rstudio-start-on-rosa.sh


### Open NEW terminal window and run:
# ssh -N -L 8000: ... (full command see instructions in first terminal window)
# Enter password


### Open http://localhost:8000 in web browser



### ===============================<<< R >>>====================================

### Set and check working directory
setwd("~/meg25/code")
getwd()


### Read in LD stats before and after thinning
before <- read_tsv("R2_before_thinning.tsv") %>%
  rename(r2 = `R^2`) %>%
  mutate(Set = "Un") %>%
  select(r2, Dprime, Set)

after <- read_tsv("R2_after_thinning.tsv") %>%
  rename(r2 = `R^2`) %>%
  mutate(Set = "Kb") %>%
  select(r2, Dprime, Set)


### Plot r2 and D statistics
boxplot(before$r2, after$r2, 
        names = c("before", "after"), 
        ylab = "r2",
        ylim = c(0, 0.1),
        outline = FALSE)   # without outliers



### =============================<<< bash >>>===================================

### Establish ssh connection to cluster in new terminal window
ssh user1234@rosa.hpc.uni-oldenburg.de
# When prompted, enter password (invisible)


### ============================================================================
### Exercise 3: Assess genetic diversity

### Calculate heterozygosity and Fis for each individual
vcftools \
  --gzvcf snps_hamlets_mac2.vcf.gz \
  --het \
  --stdout > Het_hamlets_mac2.tsv


### Measure nucleotide divergence per species / population
for i in pop_*.txt; do
vcftools \
  --gzvcf snps_hamlets_mac2.vcf.gz \
  --keep "$i" \
  --site-pi \
  --out Pi_"$i"
done


### ===============================<<< R >>>====================================

### Return to R Studio

### Read in TSV file and add population information
het <- read_tsv("Het_hamlets_mac2.tsv") %>%
  mutate(Species = str_sub(INDV, -6, -4),
         Location = str_sub(INDV, -3, -1),
         Population = str_sub(INDV, -6, -1))


### Summarize and visualize with boxplot
h <- ggplot() +
    geom_boxplot()


### Summarize and visualize with boxplot
(h <- ggplot(het, aes(x = Population, y = F, fill = Species)) +
    geom_boxplot(color = "grey20",
                 alpha = 0.75,
                 lwd = 0.3) +
    scale_fill_manual(values = c("goldenrod2", "royalblue3", "grey20", "coral2", "grey80")) +
    labs(title = NULL,
         x = "Population",
         y = "Mean genome-wide Fis") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title.y = element_text(vjust = 2),
          axis.text.x = element_text(angle = 35)
    )
)

pi 



### ============================================================================
### Solutions:

### VCF stats
bcftools stats
grep -v '#' | wc -l


### Include only sites with minor allele count of at least 2
vcftools \
    --gzvcf snps_hamlets_lg12.vcf.gz \
    --mac 2 \
    --recode \
    --stdout | bgzip > snps_hamlets_mac2.vcf.gz


### Summarize
bcftools stats snps_hamlets_mac2.vcf.gz
#> 219354 sites


### Spot check
zcat snps_hamlets_mac2.vcf.gz | tail


### Thin data so that no two sites are within 2000 bp from each other
vcftools \
    --gzvcf data/snps_hamlets_lg12.vcf.gz \
    --thin 2000 \
    --recode \
    --stdout | bgzip > snps_hamlets_2kb.vcf.gz



### ============================================================================
### Alternative code and notes

### Subset phyps2 dataset to create example file
vcftools \
    --gzvcf phyps2e_snpsfilt.vcf.gz \
    --chr LG12 \
    --keep hamlet_subset.txt \
    --recode \
    --stdout | \
    grep -v -e 'ID=Contig' -e '##GATKCommandLine=' | \
    bgzip > hamlets_LG12_snp.vcf.gz

vcftools \
    --gzvcf ~/off/phylo2/2_genotyping/out/8_geno/phyps2_all.vcf.gz \
    --chr LG12 \
    --from-bp 1 \
    --to-bp 1000000 \
    --keep hamlet_subset.txt \
    --recode \
    --stdout | \
    grep -v -e 'ID=Contig' -e '##GATKCommandLine=' | \
    bgzip > hamlets_LG12-1M_all.vcf.gz


### Summarize in ASCII
for file in R2_before_thinning.tsv R2_after_thinning.tsv; do
  echo
  echo "$file:"
  awk '$5 != "-nan" { sum += $5; n++; vals[n] = $5 }
  END {
    if (n == 0) exit
    asort(vals)
    median = (n % 2) ? vals[int((n+1)/2)] : (vals[n/2] + vals[n/2 + 1]) / 2
    print "  Mean R²: " sum/n
    print "  Median R²: " median
    print "  Count: " n
  }' "$file"
done

### Plot in ASCII
for file in R2_before_thinning.tsv R2_after_thinning.tsv; do
  echo
  echo "==> $file <=="
  awk '
  BEGIN { bin_size=0.05 }
  /^CHR/ { next }  # skip header
  {
    if ($5 != "-nan") {
      r2 = $5 + 0
      bin = int(r2 / bin_size)
      counts[bin]++
      total++
    }
  }
  END {
    max_percent = 0
    for (i in counts) {
      percent = 100 * counts[i] / total
      if (percent > max_percent) max_percent = percent
    }
    for (i = 0; i <= int(1 / bin_size); i++) {
      percent = 100 * counts[i] / total
      bar_len = int((percent / max_percent) * 50)
      bin_label = sprintf("%.2f", i * bin_size)
      printf "%5s | %s (%.1f%%)\n", bin_label, substr("##################################################", 1, bar_len), percent
    }
  }' "$file"
  echo
done