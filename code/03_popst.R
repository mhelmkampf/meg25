### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2025                             ###
### 03. Population structure                                                 ###
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



### ===============================<<< R >>>====================================

### Set and check working directory
setwd("~/meg25/code")
getwd()



### ============================================================================
### Exercise 1: Compare amount of genetic structure in three datasets

### Load packages (installed last time, otherwise use install.packages())
library(adegenet)
library(genepop)
library(pegas)


### Read in data
caribbean <- read.genepop("../data/msats/puella_caribbean.gen", ncode = 3)
timeseries <- read.genepop("../data/msats/puella_timeseries.gen", ncode = 3)
hamlets <- read.genepop("../data/msats/hamlets_caribbean.gen", ncode = 3)


### Calculate F-statistics for each locus
pegas::Fst(as.loci(caribbean))
pegas::Fst(as.loci(timeseries))
pegas::Fst(as.loci(hamlets))


### Calculate global rho_st for each dataset using genepop package
help(Fst)

# hints: use genepop::Fst(), to differentiate this function from the one in the pegas package
# rho_st was developed for microsatellites, which are allele-size based markers


### Test for genetic differentiation (G-test)
help(test_diff)



### ============================================================================
### Exercise 2: Calculate population-specific Fst

### Load packages
install.packages("hierfstat")
install.packages("tidyverse")
library(hierfstat)
library(tidyverse)


### Population-specific Fst (beta)
betas(hamlets)


### Convert output to tidyverse data frame (tibble) -- execute step by step to follow the pipeline
b <- betas(hamlets)$betaiovl %>%                   # extract Fsts from betas object
  bind_rows() %>%                                  # convert to tibble
  pivot_longer(cols = everything()) %>%            # transform data from wide to long format
  rename("Population" = name, "Fst" = value) %>%   # rename columns
  mutate(Species = str_sub(Population, 1, 3)) %>%  # add colum with species id
  arrange(desc(Fst))                               # sort data by descending Fst


### Basic plotting with ggplot
ggplot(data = b, aes(x = Population, y = Fst)) +
  geom_bar(stat = "identity")


### Make plot prettier
ggplot(data = b, 
       aes(x = reorder(Population, -Fst),       # determine basic data structure (mapping)
           y = Fst)) +
  geom_bar(stat = "identity") +                 # use barplot geom
  labs(x = "Population") +                      # change x-axis label
  theme_minimal(base_size = 16) +               # use theme "minimal" (e.g. no frame, white background)
  theme(
    panel.grid.minor = element_blank(),         # adjust theme: remove minor grid lines
    panel.grid.major.x = element_blank(),       # adjust theme: remove major grid lines intersecting x-axis
    axis.text.x = element_text(angle = 45,      # adjust theme: change angle and position of x-axis labels
                               hjust = 1, 
                               vjust = 1.25)    
  )


### Same plot, but colored by species
ggplot(data = b, 
       aes(x = reorder(Population, -Fst),
           y = Fst,
           fill = Species)) +                   # add fill mapping
  geom_bar(stat = "identity",
           color = "gray20") +                  # set bar outline color
  scale_fill_manual(values = c("#DFDF8D", "#8B4513", "#F99729",   # define fill colors for species
                               "#22198E", "#333333", "#E48175", 
                               "#B3B3B3")) +
  labs(x = "Population") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45,
                               hjust = 1, 
                               vjust = 1.25)    
  )



### ============================================================================
### Exercise 3: Visualize population structure using PCoA

### Calculate matrix of pairwise Fst
d <- genet.dist(hamlets, method = "Nei87")   # pairwise Fst following Nei (1978)
d


### Basic plotting
p <- pcoa(as.matrix(d))


### Create tidyverse data frame (tibble) with first two axes
t <- tibble(pco1 = p$vecp[, 1], 
            pco2 = p$vecp[, 2], 
            population = colnames(as.matrix(d))) %>%
  mutate(Species = str_sub(population, 1, 3))


### Plot (and store plot in variable)
g <-
  ggplot(data = t,
         aes(x = pco1, 
             y = pco2,
             color = Species)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#DFDF8D", "#8B4513", "#F99729", 
                                "#22198E", "#333333", "#E48175", "#B3B3B3")) +
  labs(x = "PCoA 1",
       y = "PCoA 2") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "transparent")
  )
g


### Add population labels
# install.packages("ggrepel")
# library(ggrepel)

# g + geom_text_repel(aes(label = population), point.padding = 10)


### How does the Fst-based PCoA compare to the Fst barplot from exercise 2?



### ============================================================================
### Exercise 1 solution

### Calculate rho_st using Genepop for each dataset
genepop::Fst("../data/msats/puella_caribbean.gen", sizes = TRUE, outputFile = "../../work/Rst_caribbean.txt")
genepop::Fst("../data/msats/puella_timeseries.gen", sizes = TRUE, outputFile = "../../work/Rst_timeseries.txt")
genepop::Fst("../data/msats/hamlets_caribbean.gen", sizes = TRUE, outputFile = "../../work/Rst_hamlets.txt")


### Test for population structure / genetic differentiation (G-test)
test_diff("../data/msats/puella_caribbean.gen", outputFile = "../../work/Diff_caribbean.gen")
test_diff("../data/msats/puella_timeseries.gen", outputFile = "../../work/Diff_timeseries.gen")
test_diff("../data/msats/hamlets_caribbean.gen", outputFile = "../../work/Diff_hamlets.gen")
