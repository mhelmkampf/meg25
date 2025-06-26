### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2025                             ###
### 09. DNA barcoding                                                        ###
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
mkdir ~/work/barcode
cd ~/work/barcode



### ============================================================================
### Exercise 1: Extract barcode sequene from Sanger reads

### Background – PCR / sequencing primer cocktail used for COI:

# VF2_t1    tgtaaaacgacggccagtCAACCAACCACAAAGACATTGGCAC
# FishF2_t1 tgtaaaacgacggccagtCGACTAATCATAAAGATATCGGCAC

# FishR2_t1 caggaaacagctatgacACTTCAGGGTGACCGAAGAATCAGAA
# FR1d_t1   caggaaacagctatgacACCTCAGGGTGTCCGAARAAYCARAA

# Amplifies approx. 650 bp fragment in 5' region of COI
# Primers are M13-tailed to allow DNA sequencing
# (Ivanova et al. 2007, Molecluar Ecology Notes, lower case = tail)


### Everyone will be working with one of the following 11 samples:
### 10, 11, 14, 15, 19, 20, 22, 26, 27, 28, 35


### From meg25/data/barcode, copy your pair of ABI trace files to the working directory
cp ~/meg25/data/barcode/seq_10F.ab1 .
cp ~/meg25/data/barcode/seq_10R.ab1 .


### View sequence (download raw file from https://github.com/mhelmkampf and upload here:)
#> https://www.gear-genomics.com/teal/


### Extract sequence from trace files
module load EMBOSS/6.6.0-foss-2023a

seqret -sequence seq_10F.ab1 -outseq seq_10F.fas
seqret -sequence seq_10R.ab1 -outseq seq_10R.fas


### Trim sequences (adjust limits according to your case)
seqkit subseq -r 65:580 seq_10F.fas > trim_10F.fas
seqkit subseq -r 50:470 seq_10R.fas > trim_10R.fas


### Reverse complement R sequence
seqkit seq -r -p trim_10R.fas > trim_10R-rc.fas


### Merge foward and reverse sequences to consensus
#>
merger \
  -asequence <forward.fas> \
  -bsequence <reverse.fas> \
  -outfile <consensus.txt> \
  -outseq <consensus.fas>



### ============================================================================
### Exercise 2: Match sequence to Barcode of Life Data Systems (BOLD)


### Visit https://boldsystems.org


### Choose "Barcode ID" | Animal Library (Public) and Rapid Species Search


### Paste your consensus sequence and submit


### Optional: Learn more about your species at https://fishbase.de



### ============================================================================
### Exercise 3: Evaluate reliability of id result


### Go to https://portal.boldsystems.org/bin and enter BIN


### Alternatively, use https://v4.boldsystems.org/index.php
### Match sequence and select "Bin page" on results page


### How confident are you in the id based on similarity score and the distribution
### of genetic distances? Is there a "barcode gap"?



### ============================================================================
### Solutions:

### Create consensus
merger \
  -asequence trim_10F.fas \
  -bsequence trim_10R-rc.fas \
  -outfile consensus_10.txt \
  -outseq consensus_10.fas


### All consensus sequences can be found at data/barcode/co1.fas


### Species IDs, similarity score and barcode gap
# 10: Oncorhynchus keta (Chum salmon / Keta-Lachs), 100%, large barcode gap
# 11: Gadus chalcogrammus (Alaska pollock / Paz. Pollack), 100%, small barcode gap
# 14: Platichthys flesus (European flounder / Flunder) or Pleuronectes sp., 100%, no barcode gap
# 15: Gadus chalcogrammus (Alaska pollock / Paz. Pollack), 100%, see above, small barcode gap
# 19: Gadus chalcogrammus (Alaska pollock / Paz. Pollack), 100%, see above, small barcode gap
# 20: Platichthys flesus (European flounder / Flunder) or Pleuronectes sp., 100%, small barcode gap
# 22: Gadus morhua (Atlantic cod / Kabeljau), 100%, small barcode gap
# 26: Gadus morhua (Atlantic cod / Kabeljau), 100%, small barcode gap
# 27:
# 28:
# 35: