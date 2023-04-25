
#!/bin/bash

#April 24, 2023

###################################################################################################################################
# DO NOT modify the next 4 lines.                                                           #######################################
REPO=/Hierarchical-tree-dissection-framework/                                               #######################################
CLUSTERING=$REPO/CLUSTERING                                                                 #######################################
## Number of Cyclospora species in population (DO NOT MODIFY THIS PART OF THE CODE).        #######################################
number_of_species_in_USA=2.                                                                 #######################################
###################################################################################################################################



###USER CAN MODIFY THE FOLLOWING LINES:

# TELL ME THE DIRECTORY WHERE YOU PASTED AND UNZIPPED THE REPO FOLDER. Make sure you modify this line!!!
repo_location=/Users/your_name/Documents/

# What is the smallest number of partitions in each population that you expect to see?
cluster_min=20

# What is the largest number of partitions in each population that you expect to see?
cluster_max=200

# Tell me the number of threads you would like to use
number_of_threads=11

## bootstraps - number of random selections of unbiased matrices (M2) do you want to make for each population.
bootstraps=1000




############  DO NOT MODIFY BELOW THIS POINT

LOC=$repo_location$CLUSTERING

matrix_folder=$repo_location/$REPO/ensemble_matrices

###########################################################################################################################################################
######                                                  ###################################################################################################
######  WRITE DIRECTORIES FOR TMP FILES & VARIABLES     ###################################################################################################
######                                                  ###################################################################################################
###########################################################################################################################################################

cd $LOC

rm -rf TMP_REP

mkdir TMP_REP
cd TMP_REP

echo $stringency > STRINGENCY
echo $cluster_min > CLUSTER_MIN
echo $cluster_max > CLUSTER_MAX
echo $matrix_folder > MATRIX_LOCATION
echo $number_of_threads > THREADS
echo $bootstraps > BOOTSTRAPS
echo $number_of_species_in_USA > SPECIES


Rscript $LOC/CLUSTER_FINDER_V2.R

cd $LOC
rm -rf TMP_REP
cd $repo_location/$REPO

echo "TREE DISSECTION COMPLETE"




