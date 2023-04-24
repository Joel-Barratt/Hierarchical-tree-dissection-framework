# Hierarchical tree dissection framework
Framework to automate the dissection of hierarchical trees into epidemiologically-meaningful partitions.

### Running this code

>Update the "Tree_dissect.sh" script between lines 17 to 32 as you see fit. We recommend that you leave the bootstrap setting at 1000, and cluster_min and cluster_max at 5 and 200 respectively. You should really only need to modify the number of threads and the repo_location variables. You will next need to place a distance matrix in csv format within the "matrix_loc" directory. We recommend testing the code on the large matrix from the associated study describing this framework: 

>**will provide link when study is published.**
>Reviewers of the paper can refer to "Supplementary file S1.xlsx"

>The distance matrix is provided within the supplementary excel spreadsheet. The distance matrix can be found in Tab B. Copy this matrix and paste it into the matrix_loc directory as a csv file. The code should automatically detect it in there so long as the file ends in csv. Only put one matrix in that folder at a time and dont forget to delete the readme.txt file within that directory when you unzip this repo.

>Next run the following code.
```bash
bash Tree_dissect.sh
```
>If all is well, then the code should run to completion (it will take a while depending on the number of threads you choose). The final cluster memberships will be printed in the "clusters_detected" directory.

# Reference
If you find this code helpful, please cite the following study: **yet to be published**.



