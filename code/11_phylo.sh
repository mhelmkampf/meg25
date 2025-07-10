### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2025                             ###
### 09. Phylogenetic inference                                               ###
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


### Create and navigate to today's working directory
mkdir ~/work/phylo
cd ~/work/phylo



### ============================================================================
### Exercise 1: Interpreting trees

# In the provided trees (see lecture slides):
# - determine what the branch lengths represent (nothing, evolutionary change, or absolute time)
# - identify a well-supported clade (monophyletic group)
# - identify a weakly supported group
# - identify a sister-group relationship
# - identify the outgroup


### ============================================================================
### Exercise 2: Infer distance-based phylogenetic tree

### Copy fasta files to working directory
cp ~/meg25/data/phylo/*.fas .


### Combine to single file
cat Rhyncho_coi_databases.fas Rhyncho_coi_samples.fas > Rhyncho_coi_combined.fas



### ============================================================================
### Exercise 3: Infer phylogenetic tree using maximum likelihood

### Align sequences using MAFFT
module load MAFFT


