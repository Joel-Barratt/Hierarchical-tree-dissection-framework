
#!/bin/bash

#UP November 19, 2020

###################################################################################################################################
# do not modify the next 2 lines.                                                           #######################################
CYCLONE=/NEW_MODULE_3_WORK/                                                                 #######################################
CLUSTERING=$CYCLONE/CLUSTERING                                                              #######################################
###################################################################################################################################

###USER MUST MODIFY THE FOLLOWING LINES:



# TELL ME THE DIRECTORY WHERE YOU PASTED THE CYCLONE FOLDER. LITERALLY WHERE YOU UNZIPPED CYCLONE AND INTEND TO RUN/INSTALL IT.
cyclone_location=/Users/joelbarratt/Documents/


# Stringency. This must be a number between 1 and 100. Usually 100 will be too high (it will not be achievable in most cases). Go for somewhere between 95 and 99
# because this will give you a sensible answer.

#stringency=97


# What number of clusters do you want to start at?
cluster_min=20


# What number of clusters do you want to finish at?
cluster_max=200


# Tell me the number of threads you would like to use
number_of_threads=11


## NUMBER OF SPECIES TO CLUTER MIGHT BE ANOTHER FLAG TO ADD? I use this in places but its not automatic.
number_of_species_in_USA=2



## bootstraps - number of random selections of unbiased matrices do you want to make?
bootstraps=1000


















#########################################################   DO NOT MODIFY BELOW THIS POINT


LOC=$cyclone_location$CLUSTERING

matrix_folder=$cyclone_location/$CYCLONE/ensemble_matrices


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
echo $number_of_species_in_USA > SPECIES # note that selecting 2 species here does not automate the process.
#code will need to be tweak (the R code) in the event that we decide there are 3 species among our genotypes..
#For example, if we start to see C. henanensis samples come our way.



#Rscript $LOC/CLUSTER_FINDER.R

Rscript $LOC/CLUSTER_FINDER_V2.R

cd $LOC

rm -rf TMP_REP

cd $cyclone_location/$CYCLONE

echo "CLUSTERING COMPLETE!"




