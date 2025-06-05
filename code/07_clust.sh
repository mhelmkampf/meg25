### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2025                             ###
### 07. Population structure / clustering                                    ###
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
### Exercise 1: Summarize and filter VCF files

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


### Convert to BED format
mkdir bed

plink --vcf hamlets_LG12_2m2k.vcf \
    --make-bed \
    --allow-extra-chr \
    --out bed/hamlets_LG12_2m2k


# Run ADMIXTURE
for k in {2..8}
do
    admixture \
    --cv -j10 \
    bed/hamlets_LG12_2m2k.bed $k > k${k}_hamlets_LG12_2m2k.out
done


# mv *.P admcv
# mv *.Q admcv


# Print CV error
for k in {2..8}
do
    grep 'CV' k${k}_hamlets_LG12_2m2k.out \
    >> CV_k2-8.out
done


# Add sample ids to proportions
for k in {1..12}
do
    paste -d " " $base/0_metadata/ids_phyps2e.txt admcv/phyps2e_m2n1l5.${k}.Q |
        sed 's/ $//g' \
        > AdmcvProp_phyps2e_m2n1l5_k${k}.tsv
done

# Plotted with R/admcv_phyps2.R


### ===============================<<< R >>>====================================

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






### Create link to today's SNP dataset
ln -s /fs/dss/home/haex1482/share/hamlets_LG12_snp.vcf.gz hamlets_LG12_snp.vcf.gz


### How many sites (SNPs) does the example dataset have?


### VCFtools is a standard software to filter and run calculations on genetic
### variation data in VCF format: https://vcftools.github.io/man_latest.html


### Use the following code block to filter and write data to a new file
### (insert filtering parameters in lieu of ...)

vcftools \
    --gzvcf hamlets_LG12_snp.vcf.gz \
    ... \
    --recode \
    --stdout | bgzip > ...


### Referring to the manual linked above, include only sites with a "Minor Allele 
### Count" (mac) of at least 2


### Re-calculate the number of sites to confirm filtering was sucessful, and
### spot check visually that there are indeed at least 2 counts of each allele



### ============================================================================
### Exercise 2: Filter by linkage disequilibrium

### Referring to the linkage decay plot, how many basepairs (bp) apart should the 
### sites be from each other to avoid linkage disequilibrium? Use vcftools --thin
### to create such a dataset


### Calculate the correlation coefficient of linkage disequilibrium in the first
### 100000 bp of both datasets (before and after thinning)
vcftools --gzvcf hamlets_LG12_snp.vcf.gz \
  --hap-r2 \
  --stdout > R2_before_thinning.tsv

vcftools --gzvcf hamlets_LG12_snp_2kb.vcf.gz \
  --hap-r2 \
  --stdout > R2_after_thinning.tsv


### Plot R2 distribution before and after thinning, in ASCII
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



### ============================================================================
### Exercise 3: Assess genetic diversity (heterozygosity and nucleotide diversity)

### Calculate heterozygosity and Fis for each individual
vcftools \
  --gzvcf hamlets_LG12_snp_mac2.vcf.gz \
  --het \
  --stdout > Het_hamlets_snp.tsv


### From here, we could average Fis per population and plot in R


### To calculate nucleotide diversity, we need all sites, not just SNPs:
ln -s /fs/dss/home/haex1482/share/hamlets_LG12-1M_all.vcf.gz hamlets_LG12-1M_all.vcf.gz


### Copy population id files to working directory
cp ../../meg25/data/gdiv/*.txt .


### Measure nucleotide divergence per species / population
for i in pop_*.txt; do
vcftools \
  --gzvcf hamlets_LG12-1M_all.vcf.gz \
  --keep "$i" \
  --site-pi \
  --stdout > Pi_"$i".tsv
done


### Average pi per population
for i in Pi_*.tsv; do
  echo $i
  awk '$3 != "-nan" { sum += $3 ; next } END { print sum / NR }' $i
  echo
done



### ============================================================================
### Solutions:

### VCF stats
bcftools stats hamlets_LG12_snp.vcf.gz | grep 'SN'
#> 1032804 sites

zcat hamlets_LG12_snp.vcf.gz | grep -c -v '#'


### Include only sites with minor allele count of at least 2
vcftools \
    --gzvcf hamlets_LG12_snp.vcf.gz \
    --mac 2 \
    --recode \
    --stdout | bgzip > hamlets_LG12_snp_mac2.vcf.gz


### Summarize
bcftools stats hamlets_LG12_snp_mac2.vcf.gz | grep 'SN'
#> 319960 sites


### Spot check
zcat hamlets_LG12_snp_mac2.vcf.gz | tail


### Thin the filtered data so that no two sites are within 2000 bp from each other
vcftools \
    --gzvcf hamlets_LG12_snp.vcf.gz \
    --thin 2000 \
    --recode \
    --stdout | bgzip > hamlets_LG12_snp_2kb.vcf.gz



### ============================================================================
### Alternative code and notes

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


zcat hamlets_LG12-1M_all.vcf.gz | /
grep -v -e 'ID=Contig' -e '##GATKCommandLine=' /
> hamlets_LG12-1M_all2.vcf.gz