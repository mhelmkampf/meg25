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


### Print combined Fasta file to screen
#>


### Align sequences using MAFFT (server version)
### https://mafft.cbrc.jp/alignment/server/

# - Copy and paste combined Fasta file into input box
# - Select Output order: Same as input, Strategy: Auto and submit



### ============================================================================
### Exercise 3: Infer phylogenetic tree using maximum likelihood

### Align sequences using MAFFT
module load MAFFT

mafft Rhyncho_coi_combined.fas > Rhyncho_coi_combined.aln


### Perform maximum likelihood + bootstrap analysis with RAxML-NG
### (https://github.com/amkozlov/raxml-ng/wiki)
module load IQ-TREE

raxml-ng \
    --all \
    --msa Rhyncho_coi_combined.aln \
    --model GTR+G \
    --bs-trees 1000


### Identify and print log file to screen
#>


### Print best ML tree with bootstrap support values to screen (Newick format)
#>


### Visualize tree on the Interactive Tree of Life website
### https://itol.embl.de

# - Upload | copy and paste tree in Newick format into "Tree text" box
# - Root tree: click on Rhina ancylostoma branch / tip | Tree structure | Re-root the tree here
# - Increase branch width: Branch option | Line style
# - Display support values: Advanced | Bootstrap / metadata > Display | Symbol or Text (range 50-100)


### How do the distance and maximum likelihood tree differ?