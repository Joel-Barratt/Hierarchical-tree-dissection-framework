
print("Now computing the cutoff distance for Cyclospora ashfordi")

########################################################################################
# ASHFORDI - big loop below BIG LOOP MULTI-THREADED:
########################################################################################

ASHFORDI_get_cutoffs_under_first_peak = mclapply(1:bootstraps, function (z) {
#ASHFORDI_get_cutoffs_under_first_peak = mclapply(1:10, function (z) {	

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

#now lets build a distance matrix from this set of unbiased specimens and decide upon a cutoff
#new_unbiased_matrix <- matrix(, nrow = length(unbiased_list_of_specimens$Seq_ID), ncol = length(unbiased_list_of_specimens$Seq_ID))
#rownames(new_unbiased_matrix) <- c(unbiased_list_of_specimens$Seq_ID)
#colnames(new_unbiased_matrix) <- c(unbiased_list_of_specimens$Seq_ID)

#for(n in rownames(new_unbiased_matrix)){
	
#	for(o in rownames(new_unbiased_matrix)){
		
#	new_unbiased_matrix[n, o] <- matrix[n, o]
		
#	}
	
#}


#d2 <- density(new_unbiased_matrix)
#plot(d2, lwd = 3.5, col = "black")


#just_colnames <-   cbind(colnames(new_unbiased_matrix), 1:length(colnames(new_unbiased_matrix)) )
#colnames(just_colnames) <- c("Seq_ID", "Number")

#species_clusters is a dataframe generated outside the loop because it takes a while to generate and we want the loop to run quickly.
species_clusters_ASHFORDI <-  merge(unbiased_list_of_specimens, species_clusters_ASHFORDI, by.x = "Seq_ID", all.x=TRUE) 

species_clusters_ASHFORDI $Cluster <- NULL
species_clusters_ASHFORDI $Cluster_name <- NULL


#now we want to generate a matrix of distances solely for inter-species distances (ASHFORDI versus ashfordi)
#lets make two separate lists for each species and then build a matrix of interspecies differences only.
species_1 <- filter(species_clusters_ASHFORDI, value == "1")
species_2 <- filter(species_clusters_ASHFORDI, value == "2")

#now lets figure out which is ASHFORDI. It will be the larger of these two lists (ASHFORDI makes up about 2 thirds of the population)

ASHFORDI_genotypes <- ""

if(length(species_1$value) < length(species_2$value)){
	
	ASHFORDI_genotypes <- species_1
	
} else {
	
	ASHFORDI_genotypes <- species_2
	
}


#now find normalised samples from within the ASHFORDI_genotypes list:







## make a normalised patristic matrix
patristic_ASHFORDI_normalized_matrix <- matrix(nrow = length(ASHFORDI_genotypes$Seq_ID), ncol = length(ASHFORDI_genotypes$Seq_ID))
rownames(patristic_ASHFORDI_normalized_matrix) <-  ASHFORDI_genotypes$Seq_ID
colnames(patristic_ASHFORDI_normalized_matrix) <- ASHFORDI_genotypes$Seq_ID

for(p in rownames(patristic_ASHFORDI_normalized_matrix)){
	
	for(q in rownames(patristic_ASHFORDI_normalized_matrix)){
		
	patristic_ASHFORDI_normalized_matrix[p, q] <- PatristicDistMatrix[p, q]
		
	}
	
}



## make a normalised raw genetic distance matrix

raw_ASHFORDI_normalized_matrix <- matrix(nrow = length(ASHFORDI_genotypes$Seq_ID), ncol = length(ASHFORDI_genotypes$Seq_ID))
rownames(raw_ASHFORDI_normalized_matrix) <-  ASHFORDI_genotypes$Seq_ID
colnames(raw_ASHFORDI_normalized_matrix) <- ASHFORDI_genotypes$Seq_ID

for(p in rownames(raw_ASHFORDI_normalized_matrix)){
	
	for(q in rownames(raw_ASHFORDI_normalized_matrix)){
		
	raw_ASHFORDI_normalized_matrix[p, q] <- matrix[p, q]
		
	}
	
}

#raw_matrix

raw_ASHFORDI_normalized_matrix <- as.dist(raw_ASHFORDI_normalized_matrix)
raw_ASHFORDI_list_norm_matrix <- dist2list(raw_ASHFORDI_normalized_matrix)
raw_ASHFORDI_list_norm_matrix_no_self2self <- filter(raw_ASHFORDI_list_norm_matrix, col != row) #remove self to self


#sorted_raw_ASHFORDI_list_norm_matrix_no_self2self <- as.data.frame(raw_ASHFORDI_list_norm_matrix_no_self2self[order(raw_ASHFORDI_list_norm_matrix_no_self2self$value),])


#patristic_matrix
patristic_ASHFORDI_normalized_matrix <- as.dist(patristic_ASHFORDI_normalized_matrix)
patristic_ASHFORDI_list_norm_matrix <- dist2list(patristic_ASHFORDI_normalized_matrix)
patristic_ASHFORDI_list_norm_matrix_no_self2self <- filter(patristic_ASHFORDI_list_norm_matrix, col != row) #remove self to self


print(z/bootstraps)

D = sort(patristic_ASHFORDI_list_norm_matrix_no_self2self$value)[rank(raw_ASHFORDI_list_norm_matrix_no_self2self$value,ties.method="random")]
pat_plus_raw_ASHFORDI <- raw_ASHFORDI_list_norm_matrix_no_self2self
pat_plus_raw_ASHFORDI$D <- D

pat_plus_raw_ASHFORDI <- as.data.frame(pat_plus_raw_ASHFORDI)

sorted_pat_plus_raw_ASHFORDI <- pat_plus_raw_ASHFORDI[order(as.numeric(pat_plus_raw_ASHFORDI$value)),]


#length(sorted_pat_plus_raw_ASHFORDI$value)

}, mc.cores= number_of_threads)



########################################################################################
## ASHFORDI -- FINAL CALCULATIONS
########################################################################################
#get_cutoffs_under_first_peak <- this should only include a set of "intra-species distances" from a normalized matrix.

ASHFORDI_keep_only_numbers <- ASHFORDI_get_cutoffs_under_first_peak

for(z in 1:length(ASHFORDI_keep_only_numbers)){
	ASHFORDI_keep_only_numbers[[z]]$col <- NULL
	ASHFORDI_keep_only_numbers[[z]]$row <- NULL
	ASHFORDI_keep_only_numbers[[z]]$value <- NULL
	ASHFORDI_keep_only_numbers[[z]] <- as.numeric(ASHFORDI_keep_only_numbers[[z]]$D)
	#ASHFORDI_keep_only_numbers[[z]] <- as.numeric(ASHFORDI_keep_only_numbers[[z]]$value)
}

ASHFORDI_average_cutoff_under_first_peak <- (Reduce("+", ASHFORDI_keep_only_numbers)/length(ASHFORDI_keep_only_numbers))
ASHFORDI_cutoff_average_matrix <- ASHFORDI_average_cutoff_under_first_peak[round((5/100)*length(ASHFORDI_average_cutoff_under_first_peak))]
ASHFORDI_golden_cluster_distance <- round(ASHFORDI_cutoff_average_matrix, digits = 3)

print("Your distance cutoff for Cyclospora ashfordi is:")
print(ASHFORDI_golden_cluster_distance)
########################################################################################




foobar <- ASHFORDI_get_cutoffs_under_first_peak

for (k in 1:length(ASHFORDI_get_cutoffs_under_first_peak)){
	
	foobar[[k]]$col <- NULL
	foobar[[k]]$row <- NULL
	foobar[[k]]$value <- as.numeric(foobar[[k]]$value)
	foobar[[k]]$D <- NULL

}

ASHFORDI_average_cutoff_foobar <- (Reduce("+", foobar)/length(foobar))
#ASHFORDI_cutoff_average_foobar_matrix <- ASHFORDI_average_cutoff_foobar$value[round((5/100)*length(ASHFORDI_average_cutoff_foobar$value))]
#ASHFORDI_golden_cluster_foobar <- round(ASHFORDI_cutoff_average_foobar_matrix, digits = 3)


ASH_dens <- density(ASHFORDI_average_cutoff_foobar$value)
plot(ASH_dens, lwd = 4, col = "black")

#ad.test(ASHFORDI_average_cutoff_foobar$value) 
#qqnorm(ASHFORDI_average_cutoff_foobar$value, main='Normal')
#qqline(ASHFORDI_average_cutoff_foobar$value)


qqnorm(ASHFORDI_average_cutoff_foobar$value, main='Normal', lwd = 1.5, col = "black", cex=1, cex.axis=1.4, cex.lab=1.4)
qqline(ASHFORDI_average_cutoff_foobar$value, lwd = 5, col = "gray60")



