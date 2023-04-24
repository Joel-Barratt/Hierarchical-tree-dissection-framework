
number_of_species = as.numeric(readLines("SPECIES"))

## The ranges below are for identifying iterations of M2. This means the tree will be dissected at 5 clusters through 1000 clusters for each population to find
## which value results in partitions where each partition has 2 or more genotypes/isolates in it.
min_cluster_range <- 5
max_cluster_range <- 1000

bootstraps = as.numeric(readLines("BOOTSTRAPS"))

setwd(matrix_location)


setwd("../CLUSTERING/TMP_REP/")
number_of_threads = as.numeric(readLines("THREADS"))



print(paste0("Now testing a range of cluster numbers from ", min_cluster_range, " to ", max_cluster_range, " to identify the maximum cluster number where no singletons exist in their own cluster."))

## Cluster distance matrix using wards method.
Ensemble_x <- as.hclust(agnes(x=matrix, diss = TRUE, stand = TRUE, method = "ward"))

pick_min_cluster_number <- c()

###just_for strongy data we can use max_cluster_range <- 200

for (k in min_cluster_range:max_cluster_range){

cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k=k))))

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

# now we need a loop to determine the number of genotypes in the smallest cluster.
# This  will usually be 2 genotypes.

cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k= number_of_clusters))))
size_of_each_cluster <- vector()

for(p in 1:number_of_clusters){

current_cluster <- filter(cluster, value %in% p)

size_of_each_cluster <- append(size_of_each_cluster,length(current_cluster$value))

}

smallest_cluster_size <- min(size_of_each_cluster)


species_clusters <- as.data.frame(melt(factor(cutree(Ensemble_x, k=number_of_species)))) # this is from the original large matrix - the Ensemble_x variable.
species_clusters <- cbind(species_clusters, c(rownames(species_clusters)))
colnames(species_clusters) <- c("value", "Seq_ID") ## this will be used later.
species_clusters_ASHFORDI <- species_clusters
species_clusters_CAYETANENSIS <- species_clusters


#extract patristic distances from your tree Ensemble_x
PatristicDistMatrix <-as.matrix(cophenetic(Ensemble_x))




