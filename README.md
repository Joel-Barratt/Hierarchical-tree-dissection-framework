# Hierarchical tree dissection framework
Framework to automate the dissection of hierarchical trees into epidemiologically-meaningful partitions.

### Running this code

>Update the "Tree_dissect.sh" script between lines 17 to 32 as you see fit. We recommend that you leave the bootstrap setting at 1000, and cluster_min and cluster_max at 5 and 200 respectively. You should really only need to modify the number of threads and the repo_location variables. You will next need to place a distance matrix in csv format within the "matrix_loc" directory. We recommend testing the code on the large matrix from the associated study describing this framework: **will provide link when study is published.**
>Reviewers of the paper can refer to "Supplementary file S1.xlsx" -- look at Tab B specifically. Copy this distance matrix and convert to CSV format.

>The distance matrix is provided within the supplementary excel spreadsheet. The distance matrix can be found in Tab B. Copy this matrix and paste it into the matrix_loc directory as a csv file. The code should automatically detect it in there so long as the file ends in csv. Only put one matrix in that folder at a time and dont forget to delete the readme.txt file within that directory when you unzip this repo.

>Next run the following code.
```bash
bash Tree_dissect.sh
```
>If all is well, then the code should run to completion (it will take a while depending on the number of threads you choose). The final cluster memberships will be printed in the "clusters_detected" directory.


## Important note
The way this code differentiates between the two species of Cyclospora is due to the fact that Cyclospora ashfordi is less common in our population than C. cayetanensis. So after the tree is cut into two, the largest of the two populations is going to be C. cayetanensis and the smaller is going to be C. ashfordi. We take this for granted as a simple way to automate the differentiation of these two types of Cyclospora, but if the proportions of each type within the population change in any way, then this part of the code will need to be edited. This step is controlled in lines 64 to 76 of the "CUTOFF_ashfordi.R" and "CUTOFF_cayetanensis.R" scripts. Also note, that one of the isolates included in our test dataset (available with the manuscript describing this method) actually belongs to the species Cyclospora henanensis (genotype C_ChHenan). The genotyping markers used as input for computation of our genetic distance matrix fail to adequately distinguish these species, and so the present code will actually result in C_ChHenan being placed within a cluster alongside isolates of C. ashfordi. 

Regarding dependencies, just open the "CLUSTER_FINDER_V2.R" script and look at the list of libraries at the top of this script to determine which packages/libraries are required to run this code.

# Reference
If you find this code helpful, please cite the following study: **yet to be published**.



