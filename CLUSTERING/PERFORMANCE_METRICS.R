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

setwd("/Users/joelbarratt/Documents/NEW_MODULE_3_WORK/clusters_detected")

number_of_threads <- 10

epi_and_clusters <- read.table("final_clusters_for_paper_PYTHON_MATRIX_repeat_OLD_METHOD_GOOD_with_epi.txt", sep="\t", header=TRUE)
#epi_and_clusters <- read.table("DYNAMIC_CUT_HYBRID.txt", sep="\t", header=TRUE)

#epi_and_clusters <- read.table("DYNAMIC_CUT_TREE_ONLY.txt", sep="\t", header=TRUE)

output_filename <- "_STRINGENCY_METHOD_REPEAT_PYTHON_MATRIX"

#output_filename <- "_NEW_METHOD_PYTHON_MATRIX"
#output_filename <- "_OLD_METHOD"
#output_filename <- "_DYNAMIC_CUT_TREE_ONLY"
#output_filename <- "_DYNAMIC_CUT_HYBRID"


## EXCLUDE ANYTHING WITHOUT AN EPI-LINK
epi_and_clusters_only_linked <- filter(epi_and_clusters, Epidemiologic.links != "Unknown epi-linkage")
unique_epidemiologic_cluster <-  unique(epi_and_clusters_only_linked $Epidemiologic.links)



### REMOVE SINGLETON EPIDEMIOLOGIC LINKS (EPI-CLUSTERS WITH ONE EPI-LINKED SAMPLE IN THEM)
Epi_clus = mclapply(unique_epidemiologic_cluster, function (o) {

epi_clusters <- filter(epi_and_clusters_only_linked, Epidemiologic.links == o)

	}, mc.cores= number_of_threads)	
	
	
#df <- bind_rows(Epi_clus_no_singles) # this will merge your list of data frames into one list instantly.
	
	
Epi_clus_no_singles	<- list()
data_frames_to_null <- c()

for(p in 1:length(Epi_clus)){
	
	if(length(rownames(Epi_clus[[p]])) != 1){
		
		Epi_clus_no_singles[[p]] <- Epi_clus[[p]]
	} else {
		
		data_frames_to_null <- append(p, data_frames_to_null)
	}
	
}


No_singleton_epi_clus <- Epi_clus_no_singles

for (z in data_frames_to_null){
	
	No_singleton_epi_clus[[z]] <- NULL
	No_singleton_epi_clus <- No_singleton_epi_clus
}


#df <- bind_rows(No_singleton_epi_clus) # this will merge your list of data frames into one list instantly.


getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

 
 for(r in 1:length(No_singleton_epi_clus)){
 	
 	#r <- 34
 	
 	No_singleton_epi_clus[[r]]$Mode <- getmode(No_singleton_epi_clus[[r]]$Cluster)	
 }



 merge_epi_with_same_clus_number <- data.frame()
 WITH_MODE_Epi_clus_no_singles <- No_singleton_epi_clus
 
 for(s in 1:length(WITH_MODE_Epi_clus_no_singles)){
	
	merge_epi_with_same_clus_number <- rbind(WITH_MODE_Epi_clus_no_singles[[s]], merge_epi_with_same_clus_number)
	
}





### REMOVE SINGLETON EPIDEMIOLOGIC LINKS (EPI-CLUSTERS WITH ONE EPI-LINKED SAMPLE IN THEM)
Epi_MODE = mclapply(unique(merge_epi_with_same_clus_number$Mode), function (t) {

MODE <- filter(merge_epi_with_same_clus_number, Mode == t)

	}, mc.cores= number_of_threads)	
	

#df <- bind_rows(Epi_MODE) # this will merge your list of data frames into one list instantly.



## NOW WE HAVE A LIST OF DATA FRAMES WITH True positives listed
for (u in 1:length(Epi_MODE)){
	
	TPs <- filter(Epi_MODE[[u]], Cluster == median(Epi_MODE[[u]]$Mode))
	
	Epi_MODE[[u]]$TP <- length(TPs$Cluster)
	
}

## NOW lets tally up the false negatives

for (w in 1:length(Epi_MODE)){
		
	Epi_MODE[[w]]$FN <- (length(rownames(Epi_MODE[[w]])) - Epi_MODE[[w]]$TP[1])
	
}

 ## NOW lets tally up the false positives.
 
# merge_epi_with_same_clus_number
 
for (x in 1:length(Epi_MODE)){

	all_in_this_cluster <- filter(merge_epi_with_same_clus_number, Cluster == median(Epi_MODE[[x]]$Mode))	
		
	Epi_MODE[[x]]$FP <- (length(rownames(all_in_this_cluster)) - Epi_MODE[[x]]$TP[1])
	
}



 ## NOW lets tally up the true negatives
 
for (y in 1:length(Epi_MODE)){

	Epi_MODE[[y]]$TN <-(length(rownames(merge_epi_with_same_clus_number)) - sum(Epi_MODE[[y]]$TP[1], Epi_MODE[[y]]$FN[1], Epi_MODE[[y]]$FP[1]))
	
}

#Now lets include epi-cluster size

for (z in 1:length(Epi_MODE)){

	Epi_MODE[[z]]$Epi_cluster_size <- length(rownames(Epi_MODE[[z]]))
	
}



## Now lets count them up


TP = mclapply(1:length(Epi_MODE), function (t) {

Epi_MODE[[t]]$TP[1]

	}, mc.cores= number_of_threads)
	
	################
	
TN = mclapply(1:length(Epi_MODE), function (t) {

Epi_MODE[[t]]$TN[1]

	}, mc.cores= number_of_threads)
	################
	
	
FP = mclapply(1:length(Epi_MODE), function (t) {

Epi_MODE[[t]]$FP[1]

	}, mc.cores= number_of_threads)
	################
		
FN = mclapply(1:length(Epi_MODE), function (t) {

Epi_MODE[[t]]$FN[1]

	}, mc.cores= number_of_threads)


write_to_file  <- bind_rows(Epi_MODE) # this will merge your list of data frames into one list instantly.

			
			
## totals		
TP <- Reduce("+", TP)
TN <- Reduce("+", TN)
FP <- Reduce("+", FP)
FN <- Reduce("+", FN)

# Now lets test performance

Sensitivity <- TP/(TP+FN)
Specificity <- TN/(TN+FP)
PPV <- TP/(TP+FP)
NPV <- TN/(TN+FN)
Accuracy <- (TP+TN)/(TP+TN+FP+FN)



### NOW LETS DO SIMPSONS INDEX
cluster_count <- max(epi_and_clusters$Cluster)

Genetic_clus = mclapply(1:cluster_count, function (a) {

genetics <- filter(epi_and_clusters, Cluster == a)

	}, mc.cores= number_of_threads)	
	

### Simpsons equation is D = 1 - ((summation n(n-1))/ N(N-1))
### n = number of genetic clusters (i.e., strains?)
### N = total number of genotypes included

for (b in 1:length(Genetic_clus)){
	
	Genetic_clus[[b]]$Numerator <- length(rownames(Genetic_clus[[b]]))*(length(rownames(Genetic_clus[[b]])) - 1)
}

Simpsons_denominator <- length(rownames(epi_and_clusters))*(length(rownames(epi_and_clusters)) - 1)


Genetic_clus_numerator = mclapply(1:length(Genetic_clus), function (c) {

numer <- Genetic_clus[[c]]$Numerator[1]

	}, mc.cores= number_of_threads)	

 final_Simpsons_numerator <- Reduce("+", Genetic_clus_numerator)
 
 
 
Simpsons_D <- (1 -(final_Simpsons_numerator/Simpsons_denominator))



#output_filename <- "_NEW_METHOD"

performance_summary <- data.frame(matrix(nrow=1,ncol=11))
colnames(performance_summary) <- c("TP", "TN", "FP", "FN", "Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "Simpsons_D", "Cluster_count")
rownames(performance_summary) <- paste0("Performance", output_filename)
performance_summary$TP <- TP
performance_summary$TN <- TN
performance_summary$FP <- FP
performance_summary$FN <- FN
performance_summary$Sensitivity <- Sensitivity
performance_summary$Specificity <- Specificity
performance_summary$PPV <- PPV
performance_summary$NPV <- NPV
performance_summary$Accuracy <- Accuracy
performance_summary$Simpsons_D <- Simpsons_D
performance_summary$Cluster_count <- max(epi_and_clusters$Cluster)


write.table(performance_summary,  paste0("PERFORMANCE_SUMMARY", output_filename,".txt"), col.names=T, quote=FALSE, sep="\t",row.names=FALSE)

write.table(write_to_file,  paste0("FINAL_PERFORMANCE_OUTPUT", output_filename,".txt"), col.names=T, quote=FALSE, sep="\t",row.names=FALSE)


setwd("/Users/joelbarratt/Documents/NEW_MODULE_3_WORK/CLUSTERING")




