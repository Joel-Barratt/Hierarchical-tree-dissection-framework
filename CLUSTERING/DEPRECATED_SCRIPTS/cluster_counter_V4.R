#UP November 19, 2020

stringency_level = as.numeric(readLines("STRINGENCY"))
number_of_threads = as.numeric(readLines("THREADS"))

date_today <- Sys.Date()


#/Users/joelbarratt/Documents///NEW_MODULE_3_WORK//ensemble_matrices
#matrix <- as.matrix(read.table("/Users/joelbarratt/Documents///NEW_MODULE_3_WORK//ensemble_matrices/2022-10-24Joel_haplotype_sheet_H_matrix.csv", sep = ",", row.names=1, header=TRUE))
## can we make this run in parallel

Ensemble <- as.matrix(matrix)
Ensemble_x <- as.hclust(agnes(x=Ensemble, diss = TRUE, stand = TRUE, metric = "manhattan", method = "ward"))


list <- dist2list(as.dist(matrix))
attach(list)
sorted_list <- list[order(value),] #### sorts by the distance in the list, smallest to largest.
sorted_list <- as.data.frame(sorted_list)
#sorted_list_filtered <- sorted_list %>% filter(value < 0.7) ###THIS WILL FILTER VALUES BELOW 0.7 BASED ON THE COLUMN CALLED "VALUE".



colnames(sorted_list) = NULL
colnames(sorted_list) <- c("Seq_ID", "row", "distance")


#CAYETANENSIS_golden_cluster_distance
#CAYETANENSIS #Dave says there should be about 28 clusters
print("Checking the number of Cyclospora cayetanensis subclusters.")

##################################################################################################################################
##################################################################################################################################

CAYETANENSIS_check_all_cluster_numbers = mclapply(min_cluster_number:max_cluster_number, function (f) {

#CAYETANENSIS_check_all_cluster_numbers = mclapply(min_cluster_number:40, function (f) {

cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k=f))))


colnames(cluster) <- NULL
cluster <- cbind(rownames(cluster), cluster)
colnames(cluster) <- c("Seq_ID", "clus_A")


cluster_first_col <-merge(sorted_list, cluster, by=c("Seq_ID"), all.x=TRUE)
colnames(cluster) <- NULL
colnames(cluster) <- c("row", "clus_B")
cluster_second_col <-merge(cluster_first_col, cluster, by=c("row"), all.x=TRUE)
cluster_calc <- subset(cluster_second_col, clus_A == clus_B)

remove_self <- cluster_calc[cluster_calc $row != cluster_calc $Seq_ID, ] # removes all self-to-self distances

distances_below_gold_cluster_distance <- remove_self %>% filter(distance < CAYETANENSIS_golden_cluster_distance)

print(f/max_cluster_number)

percent_within_cluster_distances_meeting_cutoff <- ((length(distances_below_gold_cluster_distance$distance))/length(remove_self$distance))*100 #this final result will be assigned as the output of the function.

}, mc.cores= number_of_threads)


#stringency_level <- 96
#stringency_level <- 97

names(CAYETANENSIS_check_all_cluster_numbers) <- c(paste0("clus_",min_cluster_number:max_cluster_number))


for (f in names(CAYETANENSIS_check_all_cluster_numbers)){


					if (stringency_level > CAYETANENSIS_check_all_cluster_numbers[[f]]) {
						
                                fix <- as.numeric(str_remove(f, "clus_"))
								print(paste("Searching for the most appropriate cluster number -",fix,"clusters is too small."))
	
												} else {
													
								fix <- as.numeric(str_remove(f, "clus_"))
								name_clus <- paste(date_today,"_there_are_",fix,"_genetic_clusters.txt", sep="")
								cluster2 <- as.data.frame(melt(factor(cutree(Ensemble_x, k=fix))))
								colnames(cluster2) <- NULL
								cluster2 <- cbind(rownames(cluster2), cluster2)
								colnames(cluster2) <- c("Seq_ID", "Assigned_cluster")


							write.table(cluster2, name_clus, sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
							
						print(paste("The most likely number of clusters is:", fix))
						
										if (CAYETANENSIS_check_all_cluster_numbers[[f]] >= stringency_level) {
											
											break
										}
					    }
	}


CAYETANENSIS_correct_number_of_clusters <- fix

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

ASHFORDI_check_all_cluster_numbers = mclapply(min_cluster_number:max_cluster_number, function (f) {

#CAYETANENSIS_check_all_cluster_numbers = mclapply(min_cluster_number:40, function (f) {

cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, k=f))))


colnames(cluster) <- NULL
cluster <- cbind(rownames(cluster), cluster)
colnames(cluster) <- c("Seq_ID", "clus_A")


cluster_first_col <-merge(sorted_list, cluster, by=c("Seq_ID"), all.x=TRUE)
colnames(cluster) <- NULL
colnames(cluster) <- c("row", "clus_B")
cluster_second_col <-merge(cluster_first_col, cluster, by=c("row"), all.x=TRUE)
cluster_calc <- subset(cluster_second_col, clus_A == clus_B)

remove_self <- cluster_calc[cluster_calc $row != cluster_calc $Seq_ID, ] # removes all self-to-self distances

distances_below_gold_cluster_distance <- remove_self %>% filter(distance < ASHFORDI_golden_cluster_distance)

print(f/max_cluster_number)

percent_within_cluster_distances_meeting_cutoff <- ((length(distances_below_gold_cluster_distance$distance))/length(remove_self$distance))*100 #this final result will be assigned as the output of the function.

}, mc.cores= number_of_threads)




#stringency_level <- 96
#stringency_level <- 97

names(ASHFORDI_check_all_cluster_numbers) <- c(paste0("clus_",min_cluster_number:max_cluster_number))


for (f in names(ASHFORDI_check_all_cluster_numbers)){


					if (stringency_level > ASHFORDI_check_all_cluster_numbers[[f]]) {
						
                                fix <- as.numeric(str_remove(f, "clus_"))
								print(paste("Searching for the most appropriate cluster number -",fix,"clusters is too small."))
	
												} else {
													
								fix <- as.numeric(str_remove(f, "clus_"))
								name_clus <- paste(date_today,"_there_are_",fix,"_genetic_clusters.txt", sep="")
								cluster2 <- as.data.frame(melt(factor(cutree(Ensemble_x, k=fix))))
								colnames(cluster2) <- NULL
								cluster2 <- cbind(rownames(cluster2), cluster2)
								colnames(cluster2) <- c("Seq_ID", "Assigned_cluster")


							write.table(cluster2, name_clus, sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
							
						print(paste("The most likely number of clusters is:", fix))
						
										if (ASHFORDI_check_all_cluster_numbers[[f]] >= stringency_level) {
											
											break
										}
					    }
	}


ASHFORDI_correct_number_of_clusters <- fix





#tree height to cut to get 37 clusters (as predicted using our framework) is 0.5899 and the distance threshold was 0.158
#can we use this information to determine at what height we can cut the tree for ashfordi to give the right cluster number?
#the simple multiplier or do we require an equation of a line?

## 168 clusters for ashfordi with a cutoff of 0.071 or 0.072
## 37 clusters for cayetanensis with a cutoff of about 0.158
#this is using a stringency of 97

#0.5899/0.158 = 3.733544303797468

#H <- 0.071*3.733544303797468

#cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, h=H))))

#max(as.numeric(cluster$value))

#I think equation of a line is good
#cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, h=H)))) #plot k result versus height
#determine equation of the line








#### ASHFORDI TEST COMPUTATION

#ASH_below_threshold_dist <- which(Ensemble <  0.071, arr.ind=TRUE)
#ASH_below_threshold_dist <- data.frame(ASH_below_threshold_dist)

#write.table(PatristicDistMatrix, "patristic.csv",quote = FALSE, sep = ",")
#ASH_new_df <- subset(ASH_below_threshold_dist, row != col) #FILTER OUT SELF TO SELF


#list_of_numbers_ASH <- vector(mode='list', length=length(ASH_new_df$col))
#list_of_numbers_ASH <- vector(mode='list', length=10000)

#for (k in 1:length(new_df$col)){
	
#for (k in 1:length(ASH_new_df$col)){
	
#	for (x in ASH_below_threshold_dist$row[k]){
		
#		for (y in ASH_below_threshold_dist$col[k]){
#				list_of_numbers_ASH[k] <- PatristicDistMatrix[x,y]
#				print(k/length(list_of_numbers_ASH))
#		}	
#	}
#}

#getmode <- function(ASH) {
#   uniqASH <- unique(ASH)
#   uniqASH[which.max(tabulate(match(ASH, uniqASH)))]
#}

#ASH <- unlist(list_of_numbers_ASH)
#getmode(ASH)
#median(ASH)
#mean(ASH)

#cluster_ASH <- as.data.frame(melt(factor(cutree(Ensemble_x, h= mean(ASH)/1.95))))
#max(as.numeric(cluster_ASH$value))






#### CAYETANENSIS TEST COMPUTATION

#CAY_below_threshold_dist <- which(Ensemble <  0.158, arr.ind=TRUE)
#CAY_below_threshold_dist <- data.frame(CAY_below_threshold_dist)


#CAY_new_df <- subset(CAY_below_threshold_dist, row != col)


#list_of_numbers_CAY <- vector(mode='list', length=length(CAY_new_df$col))

	
#for (k in 1:length(CAY_new_df$col)){
	
#	for (x in CAY_below_threshold_dist$row[k]){
		
#		for (y in CAY_below_threshold_dist$col[k]){
				list_of_numbers_CAY[k] <- PatristicDistMatrix[x,y]
#				print(k/length(list_of_numbers_CAY))
#		}	
#	}
#}

#getmode <- function(CAY) {
#   uniqCAY <- unique(CAY)
#   uniqCAY[which.max(tabulate(match(CAY, uniqCAY)))]
#}

#CAY <- unlist(list_of_numbers_CAY)
#getmode(CAY)
#median(CAY)
#mean(CAY)

#cluster_CAY <- as.data.frame(melt(factor(cutree(Ensemble_x, h= mean(CAY)/1.95))))
#max(as.numeric(cluster_CAY$value))











#can we use the number of clusters for cayetanensis to calculate the relative scale of our tree, so that we do not need to repeat this for ashfordi?



#cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, h=0.5899))))
#max(as.numeric(cluster$value))
#H <- (max(Ensemble_x$height)*CAYETANENSIS_golden_cluster_distance) ## need to find the multiplier
#cluster <- as.data.frame(melt(factor(cutree(Ensemble_x, h=H))))

#plot(Ensemble_x)
#abline(h=0.65,col="red")


