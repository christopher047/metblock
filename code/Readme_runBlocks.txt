### runBlocks.R  paramaters ### 

achr="chr1" ### chromosome name
input_dir="../raw/" ### directory of raw data 
min.seg=20 ### minimum number of CpG sites per segment 
bwd=300 ¤## bandwidth for distance decay function
min.block=5 ### minimum number of CpG sites per block
hclust=0.3 ### cut threshold for hclust, lower gives fewer and smaller blocks 
iqr_cutoff=10 ### (maximum-minimum) divided by the interquartile (Q3-Q1) range greater than 10 excluded 
nb=5 ### number of neighbours in impute.knn calculation
ncores=16 ### number of cores
