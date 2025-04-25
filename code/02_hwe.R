### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2025                             ###
### 02. Hardy-Weinberg equilibrium                                           ###
### ======================================================================== ###


### =============================<<< bash >>>===================================


### Establish ssh connection to cluster
ssh user1234@rosa.hpc.uni-oldenburg.de


### Create working directory for results
mkdir -p work


### Download git repository with course materials (first time only)
git clone https://github.com/mhelmkampf/meg25.git


### Update git repository (only if repository already exists)
cd meg25
git pull


### Load RStudio module
module load RStudio-Server


### Exectute script to start RStudio
rstudio-start-on-rosa.sh


### Follow instructions on screen in NEW terminal window
#> ssh -N -L 8000: ...
#> Open http://localhost:8000 in browser



### ===============================<<< R >>>====================================


### ============================================================================
### Exercise 1: Using the Genepop format


### Install and load required R packages
install.packages("adegenet")
install.packages("pegas")
library(adegenet)
library(pegas)


### Read in data from Genepop format file
help(read.genepop)
yellowblue <- read.genepop("../data/msats/yellowblue.gen", ncode = 1)


### Compare with input text file on GitHub


### Access data in the new genind object
yellowblue
yellowblue@tab
yellowblue@loc.n.all
yellowblue@all.names


### Quick test for HWE (pegas package)
hw.test(yellowblue)



### ============================================================================
### Exercise 2: Microsatellite data of Caribbean reef fish populations


### Review Genepop text file on GitHub

### Read in data from Genepop format file
barbados <- read.genepop("../data/msats/puella_barbados.gen", ncode = 3)
barbados


### What is the most / least diverse locus in terms of number of alleles?


### Locus summary using poppr (no. alleles, Simpson's index, heterozygosity, evenness)
install.packages("poppr")
install.packages("genepop")
library(poppr)
library(genepop)

locus_table(barbados)


### Is this population in HWE? What does this tell us?


### More detailed test function for HWE, with results over all alleles (genepop package)
test_HW("../data/msats/puella_barbados.gen", outputFile = "../../work/HW_barbados.txt")


### Compare these results to the Cayo de Media Luna population
medialuna <- read.genepop("../data/msats/puella_medialuna.gen", ncode = 3)
test_HW("../data/msats/puella_medialuna.gen", outputFile = "../../work/HW_medialuna.txt")


### Load and test all Caribbean populations. What patterns can you identify regarding loci and populations?
# Use data file: ../data/msats/puella_caribbean.gen



### ============================================================================
### Exercise 3 (optional): Self-fertilization and heterozygote deficit

### Hamlets are simultaneous hermaphrodites with external fertilization,
### raising the question whether self-fertilization occurs regularly

### Test for heterozygote deficiency
test_HW("../data/msats/puella_caribbean.gen", which = "global deficit", outputFile = "../../work/Hdef_caribbean.txt")



### ============================================================================
### Exercise 2 solutions

## What is the most / least diverse locus in terms of number of alleles?
barbados@loc.n.all   # Number of loci and alleles
barbados@all.names   # Names of alleles for each locus


### Is this population in HWE?
hw.test(barbados)   # Quick test


### Load and test all Caribbean populations
caribbean <- read.genepop("../data/msats/puella_caribbean.gen", ncode = 3)
test_HW("../data/msats/puella_caribbean.gen", outputFile = "../../work/HW_caribbean.txt")
