
print("Now computing the cutoff distance for Cyclospora cayetanensis")

########################################################################################
# CAYETANENSIS - big loop below - MULTI-THREADED:
########################################################################################

CAYETANENSIS_get_cutoffs_under_first_peak = mclapply(1:bootstraps, function (z) {

list_of_clusters <- list()

for (m in 1:number_of_clusters){

cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k=m))))
#generate data frame showing cluster membership.

cluster_good <- cbind(rownames(cluster), cluster)
colnames(cluster_good) <- c("Seq_ID", "Cluster")
rownames(cluster_good) <- NULL
colnames(cluster_good) <- c("Seq_ID", "Cluster")
cluster_good <- cbind(cluster_good, paste0("Cluster_",cluster_good$Cluster,"_"))
colnames(cluster_good) <- c("Seq_ID","Cluster","Cluster_name")
cluster_good <- cluster_good[order(cluster_good$Cluster),]

## generate list of data frames, where each data frame contains members of the same cluster
number_of_clusters <- max(as.numeric(cluster$value))
list_of_clusters <- list()

for(i in 1:number_of_clusters){
	
foobar <- cluster_good %>% filter(str_detect(Cluster_name, paste0("Cluster_",i,"_")))
	list_of_clusters[[paste0("Cluster_",i,"_")]] <- foobar

      }	
   }

### randomly sample two specimens (or a number representing the smallest cluster size) from each cluster
### using previously created list of data frames - to generate a normalised set.

unbiased_list_of_specimens <- NULL

for(j in 1:number_of_clusters){
current_cluster <- paste0("Cluster_",j,"_")	
current_data_frame <- list_of_clusters[[current_cluster]]
# random sampling of a number of specimens equal to the smallest cluster size, from each cluster.
random_ones <- sample_n(current_data_frame, smallest_cluster_size) 
unbiased_list_of_specimens <- rbind(unbiased_list_of_specimens, random_ones)
	
}


#species_clusters is a dataframe generated outside the loop because it takes a while to generate and we want the loop to run quickly.
species_clusters_CAYETANENSIS <-  merge(unbiased_list_of_specimens, species_clusters_CAYETANENSIS, by.x = "Seq_ID", all.x=TRUE) 

species_clusters_CAYETANENSIS $Cluster <- NULL
species_clusters_CAYETANENSIS $Cluster_name <- NULL


#now we want to generate a matrix of distances solely for inter-species distances (cayetanensis versus ashfordi)
#lets make two separate lists for each species and then build a matrix of interspecies differences only.
species_1 <- filter(species_clusters_CAYETANENSIS, value == "1")
species_2 <- filter(species_clusters_CAYETANENSIS, value == "2")

#now lets figure out which is cayetanensis. It will be the larger of these two lists (cayetanensis makes up about 2 thirds of the population)

CAYETANENSIS_genotypes <- ""

if(length(species_1$value) > length(species_2$value)){
	
	CAYETANENSIS_genotypes <- species_1
	
} else {
	
	CAYETANENSIS_genotypes <- species_2
	
}


## make a normalised patristic matrix
patristic_CAYETANENSIS_normalized_matrix <- matrix(nrow = length(CAYETANENSIS_genotypes$Seq_ID), ncol = length(CAYETANENSIS_genotypes$Seq_ID))
rownames(patristic_CAYETANENSIS_normalized_matrix) <-  CAYETANENSIS_genotypes$Seq_ID
colnames(patristic_CAYETANENSIS_normalized_matrix) <- CAYETANENSIS_genotypes$Seq_ID

for(p in rownames(patristic_CAYETANENSIS_normalized_matrix)){
	
	for(q in rownames(patristic_CAYETANENSIS_normalized_matrix)){
		
	patristic_CAYETANENSIS_normalized_matrix[p, q] <- PatristicDistMatrix[p, q]
		
	}
	
}


colnames(matrix) <- rownames(matrix)

## make a normalised raw genetic distance matrix

raw_CAYETANENSIS_normalized_matrix <- matrix(nrow = length(CAYETANENSIS_genotypes$Seq_ID), ncol = length(CAYETANENSIS_genotypes$Seq_ID))
rownames(raw_CAYETANENSIS_normalized_matrix) <-  CAYETANENSIS_genotypes$Seq_ID
colnames(raw_CAYETANENSIS_normalized_matrix) <- CAYETANENSIS_genotypes$Seq_ID

for(p in rownames(raw_CAYETANENSIS_normalized_matrix)){
	
	for(q in rownames(raw_CAYETANENSIS_normalized_matrix)){
		
	raw_CAYETANENSIS_normalized_matrix[p, q] <- matrix[p, q]
		
	}
	
}

#raw_matrix

raw_CAYETANENSIS_normalized_matrix <- as.dist(raw_CAYETANENSIS_normalized_matrix)
raw_CAYETANENSIS_list_norm_matrix <- dist2list(raw_CAYETANENSIS_normalized_matrix)
raw_CAYETANENSIS_list_norm_matrix_no_self2self <- filter(raw_CAYETANENSIS_list_norm_matrix, col != row) #remove self to self



#patristic_matrix
patristic_CAYETANENSIS_normalized_matrix <- as.dist(patristic_CAYETANENSIS_normalized_matrix)
patristic_CAYETANENSIS_list_norm_matrix <- dist2list(patristic_CAYETANENSIS_normalized_matrix)
patristic_CAYETANENSIS_list_norm_matrix_no_self2self <- filter(patristic_CAYETANENSIS_list_norm_matrix, col != row) #remove self to self


print(z/bootstraps)

D = sort(patristic_CAYETANENSIS_list_norm_matrix_no_self2self$value)[rank(raw_CAYETANENSIS_list_norm_matrix_no_self2self$value,ties.method="random")]
pat_plus_raw_CAYETANENSIS <- raw_CAYETANENSIS_list_norm_matrix_no_self2self
pat_plus_raw_CAYETANENSIS$D <- D

pat_plus_raw_CAYETANENSIS <- as.data.frame(pat_plus_raw_CAYETANENSIS)

sorted_pat_plus_raw_CAYETANENSIS <- pat_plus_raw_CAYETANENSIS[order(as.numeric(pat_plus_raw_CAYETANENSIS$value)),]



}, mc.cores= number_of_threads)



########################################################################################
## Cayetanensis -- FINAL CALCULATIONS
########################################################################################
#get_cutoffs_under_first_peak <- this should only include a set of "intra-species distances" from a normalized matrix.

CAYETANENSIS_keep_only_numbers <- CAYETANENSIS_get_cutoffs_under_first_peak

for(z in 1:length(CAYETANENSIS_keep_only_numbers)){
	CAYETANENSIS_keep_only_numbers[[z]]$col <- NULL
	CAYETANENSIS_keep_only_numbers[[z]]$row <- NULL
	CAYETANENSIS_keep_only_numbers[[z]]$value <- NULL
	CAYETANENSIS_keep_only_numbers[[z]] <- as.numeric(CAYETANENSIS_keep_only_numbers[[z]]$D)
}

CAYETANENSIS_average_cutoff_under_first_peak <- (Reduce("+", CAYETANENSIS_keep_only_numbers)/length(CAYETANENSIS_keep_only_numbers))
CAYETANENSIS_cutoff_average_matrix <- CAYETANENSIS_average_cutoff_under_first_peak[round((5/100)*length(CAYETANENSIS_average_cutoff_under_first_peak))]
CAYETANENSIS_golden_cluster_distance <- round(CAYETANENSIS_cutoff_average_matrix, digits = 3)

print("Your distance cutoff for Cyclospora cayetanensis is:")
print(CAYETANENSIS_golden_cluster_distance)
########################################################################################

