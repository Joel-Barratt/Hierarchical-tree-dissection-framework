#########LARGE DATASET ANALYSIS VERSION 2
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(dbscan)
library(gridExtra)  
library(stringr)
library(cluster)
library(phylogram)
library(msa)
library(ggtree)
library(colorspace)
library(ape)
library(seqinr)
library(colorspace)
library(phangorn)
library(randomcoloR)

library(dynamicTreeCut)


#matrix <- as.matrix(read.table("2022-10-24Joel_haplotype_sheet_H_matrix.csv", sep = ",", row.names=1, header=TRUE))
#Ensemble_x <- as.hclust(agnes(x=matrix, diss = TRUE, stand = TRUE, metric = "manhattan", method = "ward"))

Ensemble_y <- as.phylo(Ensemble_x)
tip_labels <- as.data.frame(Ensemble_y$tip.label)
colnames(tip_labels) <- "Seq_ID"
tip_labels$Tip <- rownames(tip_labels)
new_tips <- merge(cay_ash_clusters, tip_labels, by=c("Seq_ID"), all.x=TRUE)
#new_tips <- merge(dyanmic_cut, tip_labels, by=c("Seq_ID"), all.x=TRUE)



sorted_new_tips <- new_tips[order(as.numeric(new_tips$Tip)),]
sorted_new_tips$NEW_NAME <- paste0(sorted_new_tips$Seq_ID," ", sorted_new_tips$Species)
Ensemble_y$tip.label <- sorted_new_tips$NEW_NAME

cols <- c("black","gray60", "black")
Ensemble_y <- groupClade(Ensemble_y, .node=c(2114, 2113)) ###These will be the 2 clusters 2114 is ashfordi, 2113 is cayetanensis

p <- ggtree(Ensemble_y, size = 1.4, layout = "circular", aes(color = group)) + scale_color_manual(values = cols)

#p <- ggtree(Ensemble_y, size = 1.4, layout = "circular") + geom_tiplab(color = "black", size = 0.4, offset = 0.05)


#x <- p +  geom_tiplab2(color = "black", size = 0.6, offset = 0.05) #+ ####dont forget to add this plus symbol back if you want a node label.
x <- p #+  geom_tiplab(color = "black", size = 0.2, offset = 0.5) #+ ####dont forge
#x + geom_text2(aes(subset=!isTip, label=node), hjust=-0.05, size = 2, color = "red") ###

length(unique(cay_ash_clusters$Species))






palette <- distinctColorPalette(length(unique(cay_ash_clusters$Species)))
unique_genotypes <- unique(cay_ash_clusters$Species)
	
TREE_tips = mclapply(1:length(unique(cay_ash_clusters$Species)), function (n) {

these_tips <- filter(sorted_new_tips, Species == unique(cay_ash_clusters$Species)[n])
these_tips <- as.numeric(these_tips$Tip)

	}, mc.cores= number_of_threads)	
	
	
for(m in 1:length(TREE_tips)){
	x <- x + geom_hilight(node=c(TREE_tips[[m]]), fill=palette[m], extend = 0.45, alpha = 1)
}








##dynamic cuttree
	
#result_Dynamic_cut <- cutreeDynamic(Ensemble_x,  method = "tree") # -- USES ONLY THE TREE

#result_Dynamic_cut <- cutreeHybrid(Ensemble_x,  matrix) # USING HYBRID METHOD

#max(cutreeDynamic(Ensemble_x,  method = "tree"))
	
	
#dyanmic_cut <- as.data.frame(cbind(tip_labels, result_Dynamic_cut$labels))
#dyanmic_cut <- as.data.frame(cbind(tip_labels, result_Dynamic_cut))
#dyanmic_cut$Tip <- NULL	

#colnames(dyanmic_cut) <- c("Seq_ID", "Cluster")


#write.table(dyanmic_cut,  "DYNAMIC_CUT_HYBRID_TREE_ONLY.txt", col.names=T, quote=FALSE, sep="\t",row.names=FALSE)



#new_tips_DY <- merge(dyanmic_cut, tip_labels, by=c("Seq_ID"), all.x=TRUE)
#sorted_new_tips_DY <- new_tips_DY[order(as.numeric(new_tips_DY$Tip)),]
#sorted_new_tips_DY $NEW_NAME <- paste0(sorted_new_tips_DY $Seq_ID,"_Cluster_", sorted_new_tips_DY $Cluster,"_")
#Ensemble_y$tip.label <- sorted_new_tips_DY$NEW_NAME
#p <- ggtree(Ensemble_y, size = 1.4, layout = "circular") #+ scale_color_manual(values = cols)
#x <- p +  geom_tiplab(color = "black", size = 0.2, offset = 0.05) #+ ####dont forge
#length(unique(new_tips_DY$Cluster))

#palette <- distinctColorPalette(length(unique(new_tips_DY$Cluster)))
#unique_genotypes <- unique(new_tips_DY$Cluster)
	
	

#sorted_new_tips_DY$Cluster_next <- paste0("clus",sorted_new_tips_DY$Cluster,"_")

#TREE_tips = mclapply(1:26, function (n) {

#these_tips <- filter(sorted_new_tips_DY, Cluster_next == paste0("clus",n,"_"))
#these_tips <- as.numeric(these_tips$Tip)

#	}, mc.cores= number_of_threads)	
		
	
#for(m in 1:length(TREE_tips)){
#	x <- x + geom_hilight(node=c(TREE_tips[[m]]), fill=palette[m], extend = 0.45, alpha = 1)
#}







