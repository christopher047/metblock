
### covert methyl blocks files to metilene ###
### match coverage and methylation files ###

myChrSplit <- function(input_dir="../raw", pattern="*cov*")
	{
	### create list with chrs as skeys file location as value	
	covs <- list.files(path=input_dir, pattern=pattern, full.name=T) 
        stub <- unlist(lapply(strsplit(basename(covs), "_"), "[[", 1))    
        return(setNames(covs, stub)) 
	}	


matchCovMeth <- function(input_dir="../raw", output_dir="./raw/") 
	{
	### read coverage files 
	covs        <- myChrSplit(input_dir, pattern="*cov*") 	
	mets        <- myChrSplit(input_dir, pattern="*meth*") 
	### return matched coverage methylation file dataframe 
 	temp        <- cbind(data.frame(covs), data.frame(mets))
        temp$output <- 	paste0(output_dir, rownames(temp), "_metiline.tab")
	return(temp) 
	}	

getPCT <- function(cov_file, met_file) 
	{
	### creates relative meth from cov and met input files	
	require(data.table) 
	covs <- data.frame(fread(cov_file, sep="\t"), row.names=1) 
        mets <- data.frame(fread(met_file, sep="\t"), row.names=1) 
	pct  <- mets/covs
	chrs <- gsub("\\..*", "", rownames(pct))
	pos  <- gsub(".*\\.", "", rownames(pct))
	return(data.frame("chr"=chrs, "pos"=pos, pct)) 
	}


p2 <- matchCovMeth(input_dir="../raw", output_dir="./raw/")
for(i in 1:nrow(p2)) 
	{
	temp    <- getPCT("cov_file"=p2$covs[i], "met_file"=p2$mets[i])
	temp[is.na(temp)] <- "." 
       	write.table(temp, p2$output[i], row.names=FALSE, sep="\t", quote=F) 	
	}

in_files  <- list.files(path="./raw", full.name=T)
out_files <- paste0("./results/", gsub("_metiline", "_met_result", basename(in_files))) 

### run from here ###

sink("run_metiline.cmd") 
for(i in 1:length(in_files)) 
	{
	exec <- paste("metilene -a NN -b UC -m 5", in_files[i], ">", out_files[i]) 
	cat(exec, "\n") 	
	}
sink() 

#### os command ####
cmd  <- file.path(getwd(), "run_metiline.cmd") 
system(cmd) 

in_files <- list.files(path="./results", full.name=T) 
out_files <- paste0("./formatted/", gsub("_met_result", "_formatted", basename(in_files)))
for(i in 1:length(in_files)) 
	{
	fs <- file.info(in_files[i])$size[1]
	if(fs < 1){next}
	m2 <- read.table(in_files[i], sep="\t") 
	colnames(m2) <-c("chr","start","end","qvalue", "difference", "num_CpGs", "p_MWU", "p_2D KS","mean_NN", "mean_UC") 
	write.table(m2, out_files[i], quote=F, sep="\t") 
	}


