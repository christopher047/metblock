#### create chromsomome coverage and methylated counts files  ######

my_start <- Sys.time() 
#### create object and run ###
source("blockClassDefinition.R")
runChromosome <- function(achr, input_dir="../raw/", min.seg=20, bwd=300, min.block=5, hclust=0.3, iqr_cutoff=10, nb=5, ncores=16) 
	{
	test <- createMethObject() 
	test <- test@readFiles(test, input_dir=input_dir, achr=achr) 
	test <- test@createIndex(test) 
	test <- test@getSegs(test, min.seg=min.seg) 
	test <- test@imputeKNNbySegs(test, ncores=ncores, nb=nb) 
	test <- test@getBlocks(test, bwd=bwd, min.block=min.block, hclust=hclust, ncores=ncores)
	test <- test@filterIQR(test, iqr_cutoff = iqr_cutoff) 
	test <- test@calcBlockCoverage(test, min.block=min.block, ncores=ncores) 
	res  <- test@getReturnData(test)
	if(!all(!is.na(test@index))){return(c())}
        return(res) 
	}

### accessory functions###
getResultBlocks    <- function(res){blocks  <- sort(do.call("c", lapply(1:length(res), function(i){res[[i]]$blocks})))}
getResultIndices   <- function(res){blocks  <- sort(do.call("c", lapply(1:length(res), function(i){res[[i]]$index})))}
getResultMatrices  <- function(res){ m2     <- do.call("rbind", lapply(1:length(res), function(i){res[[i]]$m2}))}
getResultSegs      <- function(res){seg     <- sort(do.call("c", lapply(1:length(res), function(i){getSegRanges(res[[i]]$segs)})))}

### run for all samples  ####
common         <- paste0("chr", c(1:22, "X", "Y", "M")) 
tot_70         <- lapply(1:length(common), function(i){runChromosome(common[i])})
names(tot_70)  <- common

cat("finished\n") 

### get results ####
blocks <- getResultBlocks(tot_70)
index  <- getResultIndices(tot_70) 
mat    <- getResultMatrices(tot_70) 
segs   <- getResultSegs(tot_70) 

cat("Saving data\n") 
#### write results ####
temp_block           <- data.frame(blocks) 
rownames(temp_block) <- names(blocks) 
write.table(temp_block, "../results/raw_blocks_trim.tab", sep="\t") 

#### write index ####
index           <- data.frame(index)
rownames(index) <- paste0(index$seqnames, ".", index$start)
write.table(index, "../results/raw_index_trim.tab", sep="\t") 

### write imputed matrix ###
write.table(mat, "../results/raw_mat_trim.tab", sep="\t") 

### write segments #####
temp <- data.frame(apply(data.frame(segs), 2, as.character))  
write.table(temp, "../results/segs_trim.tab", sep="\t") 


cat("Total time\n") 
cat(Sys.time()-my_start) 
cat("\n") 


