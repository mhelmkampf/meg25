### ======================================================================== ###
### Exercises in Marine Ecological Genetics 2025                             ###
### 04. Genome sequencing and assembly                                       ###
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



### ============================================================================
### Exercise 1: Handle genome assemblies on the command line

### Reorganize work directory
cd ~/work

mkdir msats
mv *.* msats

mkdir asm
cd asm


### Create link to hamlet genome assemblies
ln -s /fs/dss/home/haex1482/share/HypPue1_assembly_hybrid.fas HypPue1_assembly_hybrid.fas
ln -s /fs/dss/home/haex1482/share/HypPue2_assembly_pacbio.fas HypPue2_assembly_pacbio.fas


### Print assembly / first 10 lines to screen
cat HypPue1_assembly_hybrid.fas
head HypPue1_assembly_hybrid.fas


### Print only first line (header of first sequence)


### Look at a slice of underlying read data




### ============================================================================
### Exercise 2: Assess read quality before and after trimming



### ============================================================================
### Solutions

head -n 1 HypPue1_assembly_hybrid.fas


