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
### Exercise 1: 

### From meg25/data/barcode, copy your pair of ABI trace files to the working directory 
cp ~/meg25/data/barcode/seq_10F .
cp ~/meg25/data/barcode/seq_10R .





### ============================================================================
### Solutions:

