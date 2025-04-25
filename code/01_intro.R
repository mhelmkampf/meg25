### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2025                             ###
### 01. Introduction: R basics                                               ###
### ======================================================================== ###


### =============================<<< bash >>>===================================


### Establish ssh connection to cluster
ssh user1234@rosa.hpc.uni-oldenburg.de


### Download git repository with course materials
git clone https://github.com/mhelmkampf/meg25.git


### Load RStudio module
module load RStudio-Server


### Exectute script to start RStudio
rstudio-start-on-rosa.sh


### Follow instructions on screen in new terminal window
#> ssh -N -L 8000: ...
#> Open http://localhost:8000 in browser



### ===============================<<< R >>>====================================


### ============================================================================
### Exercise 1: Manual test for HWE


### What are the allele frequencies?

### What are the expected genotype frequencies?

### Is the population in Hardy-Weinberg equilibrium?


### ----------------------------------------------------------------------------
### Genotype-phenotype relationship: yy = yellow, yb = green, bb = blue phenotype


### Allele frequencies
n <- 10                        # number of individuals
y <- ((6 * 2) + 1) / (n * 2)   # yellow allele frequency
b <- ((3 * 2) + 1) / (n * 2)   # blue allele frequency


### Expected genotype frequencies
yy_e <- y ^ 2
yb_e <- 2 * y * b
bb_e <- b ^ 2


### Enter data into data frame for test
dat <- data.frame(row.names = c("yy", "yb", "bb"),
                  "observed" = c(6, 1, 3),
                  "expected" = c(yy_e, yb_e, bb_e)
)


### Perform Pearson's chi-squared test of goodness of fit
help(chisq.test)

test <- chisq.test(dat$observed, p = dat$expected)
test
pchisq(test$statistic, df = 1, lower.tail = FALSE)   # recalculate p with 1 degree of freedom

# Note: The Chi-square test is not actually recommended for such low sample sizes, 
# we use it here because it is appropriate for most real world data


### Heterozygosity
Ho <- 1 / n
He <- 2 * y * b


### Fixation index
Fis = (He - Ho) / He



### ============================================================================
### Links and resources

# R cheatsheet: https://cran.r-project.org/doc/contrib/Short-refcard.pdf
# R for Beginners: YaRrr! The Pirate's Guide to R, https://bookdown.org/ndphillips/YaRrr/
# Advanced R: R for Data Science, https://r4ds.had.co.nz/index.html
