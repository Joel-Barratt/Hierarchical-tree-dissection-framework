

## select the range of clusters to check for the maximum number of clusters that has at least 2 specimens in it (no less).
min_cluster_range <- 5
max_cluster_range <- 1000

bootstraps <- 1000 # number of times to run the code in order to calculate an average cutoff

### first we need to select and import the most recent matrix
matrix_location <- readLines("MATRIX_LOCATION")

setwd(matrix_location)
myFiles <- list.files(pattern="*csv", all.files = FALSE)
myMatrix <- (sort(myFiles, decreasing = TRUE)[1])
matrix <- as.matrix(read.table(myMatrix, sep = ",", row.names=1, header=TRUE))


setwd("../CLUSTERING/TMP_REP/")
number_of_threads = as.numeric(readLines("THREADS"))



print(paste0("Now testing a range of cluster numbers from ", min_cluster_range, " to ", max_cluster_range, " to identify the maximum cluster number where no singletons exist in their own cluster."))
## Cluster_distance_matrix
Ensemble_x <- as.hclust(agnes(x=matrix, diss = TRUE, stand = TRUE, metric = "manhattan", method = "ward"))

pick_min_cluster_number <- c()

for (k in min_cluster_range:max_cluster_range){
#for (k in 50){

cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k=k))))

#cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k=48))))
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
	
	if(length(foobar$Seq_ID) == 1){
			
		pick_min_cluster_number <- cbind(pick_min_cluster_number, print(k-1))
		
		break
	} else {
	
	list_of_clusters[[paste0("Cluster_",i,"_")]] <- foobar
	
### Now we need an "if" loop that checks if any given cluster has only one specimens. if a cluster has only 1 specimens the loop will re-cluster
### the matrix at the previous cluster number, and then it will save the list of clusters variable for the following steps.
      }	
   }
}

number_of_clusters <- min(pick_min_cluster_number[,1])



# test if this cluster number gives you a minimum of two specimens in the smallest cluster
#cluster_fooo <- as.data.frame(melt(factor(cutree(Ensemble_x, k= number_of_clusters))))
#colnames(cluster_fooo) <- NULL
#cluster_fooo <- cbind(rownames(cluster_fooo), cluster_fooo)
#colnames(cluster_fooo) <- c("Seq_ID", "Cluster")
#cluster_fooo <- cbind(cluster_fooo, paste0("Cluster_",cluster_fooo$Cluster,"_"))
#colnames(cluster_fooo) <- c("Seq_ID", "Cluster","Cluster_name")

#floop <- c()
#for(q in 1:length(unique(cluster_fooo$Cluster))){
#	temp <- cluster_fooo %>% filter(str_detect(Cluster_name, paste0("Cluster_",q,"_")) == TRUE)
#	floop <- cbind(floop, print(length(temp[,1])))
#}


#print("If the following statements are true, then things are looking good: -#")
#print("The following value should be 2 or greater:")
#print(min(floop))


# test if this cluster number gives you a minimum of two specimens in the smallest cluster
#cluster_fooo <- as.data.frame(melt(factor(cutree(Ensemble_x, k= (number_of_clusters + 1)))))
#colnames(cluster_fooo) <- NULL
#cluster_fooo <- cbind(rownames(cluster_fooo), cluster_fooo)
#colnames(cluster_fooo) <- c("Seq_ID", "Cluster")
#cluster_fooo <- cbind(cluster_fooo, paste0("Cluster_",cluster_fooo$Cluster,"_"))
#colnames(cluster_fooo) <- c("Seq_ID", "Cluster","Cluster_name")

#floopy <- c()
#for(q in 1:length(unique(cluster_fooo$Cluster))){
#	temp <- cluster_fooo %>% filter(str_detect(Cluster_name, paste0("Cluster_",q,"_")) == TRUE)
#	floopy <- cbind(floopy, print(length(temp[,1])))
#}

#print("The following value should be equal to 1:")
#print(min(floopy))
#print("If these statements are not true, then there may be a problem. If the second value is not equal to one, then a higher number of clusters should be selected to produce your unbiased dataset.")




# now we need a loop to determine the number of specimens in the smallest cluster
# just in case we can select more than 2 specimens in the next step. If we can, then we should for the sake of diversity.
# This value will usually be 2 however.

cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k= number_of_clusters))))
size_of_each_cluster <- vector()

for(p in 1:number_of_clusters){

current_cluster <- filter(cluster, value %in% p)

size_of_each_cluster <- append(size_of_each_cluster,length(current_cluster$value))

}

smallest_cluster_size <- min(size_of_each_cluster)







#### THE FOLLOWING PART YOU REPEAT A THOUSAND TIMES SO THAT YOU CAN GET AN ACCURATE AVERAGE CUTOFF FROM MANY REPLICATES
#### otherwise you may get bouncy specimens again and the epis dont really like that. Paralellize it so it doesn't take a decade.

get_average_of_cutoffs = mclapply(1:bootstraps, function (z) {

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

random_ones <- sample_n(current_data_frame, smallest_cluster_size) # random sampling of a number of specimens equal to the smallest cluster size, from each cluster.


unbiased_list_of_specimens <- rbind(unbiased_list_of_specimens, random_ones)
	
}


#now lets build a distance matrix from this set of unbiased specimens and decide upon a cutoff

new_unbiased_matrix <- matrix(, nrow = length(unbiased_list_of_specimens$Seq_ID), ncol = length(unbiased_list_of_specimens$Seq_ID))

rownames(new_unbiased_matrix) <- c(unbiased_list_of_specimens$Seq_ID)
colnames(new_unbiased_matrix) <- c(unbiased_list_of_specimens$Seq_ID)


for(n in rownames(new_unbiased_matrix)){
	
	for(o in rownames(new_unbiased_matrix)){
		
	new_unbiased_matrix[n, o] <- matrix[n, o]
		
	}
	
}





#library(cluster)
#library(phylogram)
#library(ggtree)
#unbiased <- as.matrix(new_unbiased_matrix)
#unbiased_norm <- (unbiased - min(unbiased))
#unbiased_norm <- unbiased_norm/max(unbiased_norm)
#d2 <- density(unbiased_norm)
#plot(d2, lwd = 3.5, col = "black")
#Ensemble_x <- as.phylo(as.hclust(agnes(x= unbiased_norm, diss = TRUE, stand = TRUE, metric = "manhattan", method = "ward")))
#p <- ggtree(Ensemble_x, size = 0.7, layout = "circular")
#p


list <- dist2list(as.dist(new_unbiased_matrix))

attach(list)
sorted_list <- list[order(value),] #### sorts by the distance in the list, smallest to largest.

got_your_tail <- sorted_list[sorted_list$col != sorted_list$row, ] # removes all self-to-self distances

#I ran these lines to examine what the tree and the density plot of the new unbiased matrix look like
# examination of these can be helpful to understand how this is all working
#d2 <- density(got_your_tail$value)
#plot(d2, lwd = 3.5, col = "black") + abline(v=0.19 , col="red", lty="dashed", lwd=2)

cutoff <- got_your_tail[round((5/100)*length(got_your_tail$value)),3]

}, mc.cores= number_of_threads)





### UP TO HERE. NOTE THAT THE LINE BELOW




#source("../get_all_matrices.R")

#get_average_of_cutoffs


calculate_average <- list()

for (i in 1:length(get_average_of_cutoffs)){
	
	calculate_average <- append(calculate_average, get_average_of_cutoffs[[i]])	
}


cutoff <- sum(as.numeric(calculate_average))/length(calculate_average)

cutoff <- round(cutoff, digits = 2)

print("Your distance cutoff is:")
print(cutoff)

golden_cluster_distance <- cutoff
#golden_cluster_distance_OLD_WAY <- cutoff
