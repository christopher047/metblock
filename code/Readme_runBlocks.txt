### runBlocks.R  paramaters ### 

achr="chr1" ### chromosome name
input_dir="../raw/" ### directory of raw data 
min.seg=20 ### minimum number of CpG sites per segment 
bwd=300 ¤## bandwidth for distance decay function
min.block=5 ### minimum number of CpG sites per block
hclust=0.3 ### cut threshold for hclust 
iqr_cutoff=10 ### Block range divided by block IQR threshold 
nb=5 ### number of neighbours in impute.knn calculation
ncores=16 ### number of cores
