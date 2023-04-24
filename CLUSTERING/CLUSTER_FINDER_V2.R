library(cluster)
library(tidyverse)
library(quantmod)
library(dplyr)
library(stringr)
library(filesstrings)
library(reshape)
library(dbscan)
library(spaa)
require(parallel)
require(stringr)
library(cluster)
library(tidyverse)
library(ape)
library(utils)
library(colorspace)
library(DescTools)

#UP November 19, 2020

date_today <- Sys.Date()

temporary_directory <- getwd()




#### First you need to import the most recent ensemble matrix as a variable in R called "matrix" in R

matrix_location <- readLines("MATRIX_LOCATION")

setwd(matrix_location)

myFiles <- list.files(pattern="*csv", all.files = FALSE)
myMatrix <- (sort(myFiles, decreasing = TRUE)[1]) 
matrix <- as.matrix(read.table(myMatrix, sep = ",", row.names=1, header=TRUE))


setwd(temporary_directory)

##### Now you need to calculate the number of genetic clusters using the "gold standard" list. Use this to calculate the number of genetic clusters. "number_of_genetic_clusters" needs to be the name of this variable.


###set range of cluster numbers to test
min_cluster_number = as.numeric(readLines("CLUSTER_MIN"))
max_cluster_number = as.numeric(readLines("CLUSTER_MAX"))


### below here is where things change completely to the original:


source("../INTRASPECIES_cutoffs_V1.R")


source("../CUTOFF_cayetanensis.R") #0.786
source("../CUTOFF_ashfordi.R") #0.406 or 0.407


#CAYETANENSIS_golden_cluster_distance
#ASHFORDI_golden_cluster_distance

#species_clusters
number_1 <- filter(species_clusters, value == 1)
number_2 <- filter(species_clusters, value == 2)

number_1$Species <- ""
number_2$Species <- ""

if(length(number_1$value) > length(number_2$value)){
	
	number_1$Species <- "cayetanensis"
	number_2$Species <- "ashfordi"
	
} else {
	
	number_2$Species <- "cayetanensis"
	number_1$Species <- "ashfordi"
	
}



ASHFORDI_clusters <- as.data.frame(melt(factor(cutree(Ensemble_x, h= ASHFORDI_golden_cluster_distance))))
CAYETANENSIS_clusters <- as.data.frame(melt(factor(cutree(Ensemble_x, h=CAYETANENSIS_golden_cluster_distance))))




#ASHFORDI_clusters <- as.data.frame(melt(factor(cutree(Ensemble_x, k= 84))))
#CAYETANENSIS_clusters <- as.data.frame(melt(factor(cutree(Ensemble_x, k=36))))




#p <- ggtree(Ensemble_x, size = 1.4, layout = "circular")

CAYETANENSIS_clusters$Seq_ID <- rownames(CAYETANENSIS_clusters)
ASHFORDI_clusters$Seq_ID <- rownames(ASHFORDI_clusters)


number_1$value <- NULL
number_2$value <- NULL



if(unique(number_1$Species) == "cayetanensis"){
	
	CAYETANENSIS_clusters <- merge(number_1, CAYETANENSIS_clusters, by=c("Seq_ID"), all.x=TRUE)
	ASHFORDI_clusters <- merge(number_2, ASHFORDI_clusters, by=c("Seq_ID"), all.x=TRUE)
	
} else {
	
	CAYETANENSIS_clusters <- merge(number_2, CAYETANENSIS_clusters, by=c("Seq_ID"), all.x=TRUE)
	ASHFORDI_clusters <- merge(number_1, ASHFORDI_clusters, by=c("Seq_ID"), all.x=TRUE)

}




CAYETANENSIS_clusters$value <- as.numeric(CAYETANENSIS_clusters$value)
ASHFORDI_clusters$value <- as.numeric(ASHFORDI_clusters$value)

sorted_CAYETANENSIS_clusters <- CAYETANENSIS_clusters[order(CAYETANENSIS_clusters$value),]
sorted_ASHFORDI_clusters <- ASHFORDI_clusters[order(ASHFORDI_clusters$value),]

sorted_CAYETANENSIS_clusters$Species <- paste0("cayetanensis_", sorted_CAYETANENSIS_clusters $value,"_")
sorted_ASHFORDI_clusters$Species <- paste0("ashfordi_", sorted_ASHFORDI_clusters $value, "_")
cay_ash_clusters <- rbind(sorted_CAYETANENSIS_clusters, sorted_ASHFORDI_clusters)

cay_ash_clusters$value <- NULL

#length(unique(cay_ash_clusters$Species))




source("../TREE.R")

x # will print the tree




CLUS_RENAME = mclapply(1:length(unique(cay_ash_clusters$Species)), function (s) {
these_tips <- filter(cay_ash_clusters, Species == unique(cay_ash_clusters$Species)[s])
	}, mc.cores= number_of_threads)
	
final_clusters <- data.frame()
for(t in 1:length(CLUS_RENAME)){		
	CLUS_RENAME[[t]]$Cluster <- t
	final_clusters <- rbind(CLUS_RENAME[[t]], final_clusters)	
	}	



for(u in 1:length(rownames(final_clusters))){
	
	if(grepl("cayetanensis", final_clusters[u,2]) == TRUE){
		
		final_clusters[u,2] <- "Cyclospora cayetanensis"
		
	} else {
		
		final_clusters[u,2] <- "Cyclospora ashfordi"

	}
}



#final_clusters$Species <- NULL
sorted_final_clusters <- final_clusters[order(final_clusters$Cluster),]

setwd("..")
setwd("../clusters_detected")
write.table(sorted_final_clusters, "final_clusters_for_paper_PYTHON_MATRIX.txt", sep="\t", quote=FALSE, row.names=FALSE)

#write.table(sorted_final_clusters, "final_clusters_for_paper_PYTHON_MATRIX_repeat_OLD_METHOD_GOOD.txt", sep="\t", quote=FALSE, row.names=FALSE)




####




