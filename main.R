### requirements ###
#require(impute)
#require(parallel)
#require(data.table)
#require(GenomicRanges)
#require(bsseq)
#require(BiocParallel)
#require(dmrseq) 
#require(ReactomePA)
#require(Gviz)
#require(rtracklayer)
#require(TxDb.Hsapiens.UCSC.hg38.knownGene)
#require(org.Hs.eg.db)
#require(clusterProfiler)
#requirt(R.utils) 

#metilene must be in path 
#raw data used for all three methods is in folder raw

start_time <- Sys.time()
base       <- getwd()
set.seed(1024) 

### unpack data ##
setwd(file.path(base, "code"))
source("unpack_data.R")
setwd(base)

### run blocks ###
setwd(file.path(base, "code")) 
source("runBlocks.R") 
setwd(base) 

### metilene ###
setwd(file.path(base, "metilene")) 
source("mb3met.R") 
setwd(base) 

### dmrseq  ###
setwd(file.path(base, "dmrseq"))
source("RunDMRSeq.R") 
setwd(base) 

### tables ###
setwd(file.path(base, "code")) 
source("masterTable.R") 
setwd(base) 

### rename to match publication ###
setwd(file.path(base, "code"))
source("renameTables.R")
setwd(base)

end_time <- Sys.time()-start_time
cat(paste("finished", end_time, "\n")) 

