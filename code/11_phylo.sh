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
cat 




### ============================================================================
### Exercise 2: Match sequence to Barcode of Life Data Systems (BOLD)


### Visit https://v4.boldsystems.org/index.php
### Note: https://boldsystems.org (version 5) does not have full functionality yet


### Select "Identification v4" | Animal Identification (COI) | Species Level Barcode Records (default)


### Paste your consensus sequence and submit


### Explore the results:
### - What species has been identified?
### - Are there multiple species in the top 20? What is the similarity distribution?
### – Is it a plausible match considering the samples were obtained from seafood?


### Optional: Learn more about your species at https://fishbase.de


### Evaluate the reliability of your result by selecting its "Bin page"
### - What is the maximum p-distance within the species?
### – How does the distance distribution compare to the nearest neighbor? Is there a "barcode gap"?
### - How confident are you overall in your id?


### Optional: Identify barcode with BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi



### ============================================================================
### Solutions:

### Create consensus
merger \
  -asequence trim_10F.fas \
  -bsequence trim_10R-rc.fas \
  -outfile consensus_10.txt \
  -outseq consensus_10.fas


### All consensus sequences can be found in data/barcode/co1.fas


### Species IDs, BIN, barcode gap

# 10: Oncorhynchus keta (Chum salmon / Keta-Lachs), BOLD:AAA3872, large barcode gap
# 11: Gadus chalcogrammus (Alaska pollock / Paz. Pollack), BOLD:AAA3069, small barcode gap
# 14: Platichthys flesus (European flounder / Flunder) or Pleuronectes sp., cannot be identified to species
# 15: Gadus chalcogrammus (Alaska pollock / Paz. Pollack), BOLD:AAA3069, small barcode gap
# 19: Gadus chalcogrammus (Alaska pollock / Paz. Pollack), BOLD:AAA3069, small barcode gap
# 20: Platichthys flesus (European flounder / Flunder) or Pleuronectes sp., cannot be identified to species
# 22: Gadus morhua (Atlantic cod / Kabeljau) or Gadus sp., small barcode gap / cannot be identified to species
# 26: Gadus morhua (Atlantic cod / Kabeljau) or Gadus sp. small barcode gap /cannot be identified to species
# 27: Pollachius virens (Pollock / Köhler), BOLD:AAB2980, large barcode gap
# 28: Oncorhynchus keta (Chum salmon / Keta-Lachs), BOLD:AAA3872, large barcode gap
# 35: Coregonus sp. (Cisco), cannot be identified to species