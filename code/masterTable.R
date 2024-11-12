### load annotation ###
source("getGTF.R") 
ucsc <- getGTF()
reg  <- getRegSites()

### et qvalue for metiline and dmrseq
qvalue <- 0.05

### get block processed data ###
#raw_blocks_T70_H4.tab  raw_index_T70_H4.tab  raw_mat_T70_H4.tab  segs_T70_H4.tab
mat    <- read.table("../results/raw_mat_trim.tab", header=T, row.names=1, sep="\t")
index  <- GRanges(read.table("../results/raw_index_trim.tab", header=T, row.names=1, sep="\t")) 
blocks <- GRanges(read.table("../results/raw_blocks_trim.tab", header=T, row.names=1, sep="\t")) 

### get metilene data ####
met_files <- list.files(path="../metilene/formatted/", pattern="formatted.tab", full.name=T)
met_files <- met_files[which(file.size(met_files)>0)]
temp      <- lapply(1:length(met_files), function(i){read.table(met_files[i], header=T, sep="\t")})
mets      <- GRanges(data.frame(do.call("rbind", temp))) 
mets      <- sort(mets) 
sig_mets  <- mets[mets$qvalue <= qvalue]

### get dmrseq data ###
dmr_files <- list.files("../dmrseq/results", pattern="_regions.tab", full.name=T)			
dmr_files <- dmr_files[which(file.size(dmr_files)>0)]
temp      <- lapply(1:length(dmr_files), function(i){read.table(dmr_files[i], header=T, sep="\t")})
dmrs      <- GRanges(data.frame(do.call("rbind", temp)))
dmrs      <- sort(dmrs)
sig_dmrs  <- dmrs[dmrs$qval <= qvalue]

### get total sites in all chromosomes ### 
all_files <- list.files(path="../raw", pattern="*_cov.tab", full.name=T) 
tot_sites <- lapply(1:length(all_files), function(i){nrow(read.table(all_files[i], header=T, sep="\t"))})
tot_sites <- sum(unlist(tot_sites))

### table basic sites 1###
r1             <- c(length(blocks), length(sig_dmrs), length(sig_mets)) 
### percentage of total sites in result
r2             <- c(sum(blocks$L)/tot_sites, sum(sig_dmrs$L)/tot_sites, sum(sig_mets$num_CpGs)/tot_sites)*100   
### average block size 
r3             <- c( mean(width(blocks)), mean(width(sig_dmrs)), mean(width(sig_mets)))
### average number of CpGs in block 
r4             <- c(mean(blocks$L), mean(sig_dmrs$L), mean(sig_mets$num_CpGs))
tab1           <- do.call("rbind", list(r1, r2, r3, r4))
colnames(tab1) <- c("blocks", "dmrseq", "metilene") 
rownames(tab1) <- c("number", "pct of total", "mean width", "mean number of CpGs")
tab1           <- round(tab1, 2) 
write.table(tab1, "../tables/table1.tab", sep="\t") 

### table 2 overlaps 
tabs <- c("blocks", "dmrseq", "metline")
tab2 <- matrix(nrow=3, ncol=3, 0, dimnames=list(tabs, tabs)) 
tab2[1,2] <- sum(countOverlaps(blocks, sig_dmrs))/length(blocks) 
tab2[1,3] <- sum(countOverlaps(blocks, sig_mets))/length(blocks)
tab2[2,1] <- sum(countOverlaps(sig_dmrs, blocks))/length(sig_dmrs) 
tab2[2,3] <- sum(countOverlaps(sig_dmrs, sig_mets))/length(sig_dmrs)
tab2[3,1] <- sum(countOverlaps(sig_mets, blocks))/length(sig_dmrs)
tab2[3,2] <- sum(countOverlaps(sig_mets, sig_dmrs))/length(sig_mets) 
tab2      <- round(tab2*100,1)
tab2[1,1] <- "x"
tab2[2,2] <- "x"
tab2[3,3] <- "x"
write.table(tab2, "../tables/table2.tab", sep="\t")


### table 3 functional annotion ###
getAnnot <- function(obj,val) 
	{
	agtf        <- obj[val][[1]]		
	temp_blocks <- sum(countOverlaps(blocks, agtf)>0)/length(blocks)  
	temp_dmrs   <- sum(countOverlaps(sig_dmrs, agtf)>0)/length(sig_dmrs)  
	temp_mets   <- sum(countOverlaps(sig_mets, agtf)>0)/length(sig_mets)  
	return(list("blocks"=temp_blocks, "dmrseq"=temp_dmrs, "metiline"=temp_mets)) 
       	}	
temp           <- lapply(names(ucsc), function(n){getAnnot(ucsc, n)}) 
names(temp)    <- names(ucsc) 
part1          <- do.call("rbind", temp) 
temp           <- lapply(names(reg), function(n){getAnnot(reg,n)}) 
names(temp)    <- names(reg) 
part2          <- do.call("rbind", temp) 
temp           <- rbind(part1, part2) 
tab3           <- apply(temp, 2, as.numeric) 
rownames(tab3) <- rownames(temp) 
tab3           <- round(tab3*100,2) 
tab3           <- tab3[-c(1,2,13,14),] 
write.table(tab3, "../tables/table3.tab", sep="\t")


### table 4 gene enrichment  ####
library(ReactomePA)
all_genes <- ucsc$full[ucsc$full$type == "gene"]
fmet      <- na.omit(all_genes[queryHits(findOverlaps(all_genes, sig_mets))]$entrez) 
fmet_res  <- data.frame(enrichPathway(fmet, readable=T))
fdmrs     <- na.omit(all_genes[queryHits(findOverlaps(all_genes, sig_dmrs))]$entrez)
fdmr_res  <- data.frame(enrichPathway(fdmrs, readable=T))
fblk      <- na.omit(all_genes[queryHits(findOverlaps(all_genes, blocks))]$entrez)
fblk_res  <- data.frame(enrichPathway(fblk, readable=T))
### how to save ###
write.table(fmet_res, "../tables/metiline_genrich.tab", sep="\t") 
write.table(fdmr_res, "../tables/dmrseq_genrich.tab", sep="\t") 
write.table(fblk_res, "../tables/blocks_genrich.tab", sep="\t") 

#### table 5 only in blocks  ####
keep <- setdiff(rownames(fblk_res), union(rownames(fmet_res), rownames(fdmr_res))) 
write.table(fblk_res[keep,], "../tables/only_in_blocks.tab", sep="\t") 

### start drawing examples of multiple overlaps 
### need to rename blocks for indexing pathway results ####

ov    <- findOverlaps(blocks, sig_dmrs) 
temp  <- aggregate(queryHits(ov), list(subjectHits(ov)), c) 
temp  <- temp[order(unlist(lapply(temp$x, length)), decreasing=T),] 



#### Gviz #####
library(Gviz) 

### load meth pct from raw_files###
getRawPct <- function(achr, path="../raw/")
	{
	cov <- read.table(paste0(path, achr, "_tot_cov.tab"), header=T, row.names=1, sep="\t") 
	met <-  read.table(paste0(path, achr, "_tot_meth.tab"), header=T, row.names=1, sep="\t")
	chromosome <- gsub("\\..*", "", rownames(cov))
	end        <- gsub(".*\\.", "", rownames(cov))
	temp       <- GRanges(data.frame("chromosome"=chromosome, "start"=end, "end"=end, met/cov))  
	return(temp)  
	}

all_chrs <- unlist(lapply(strsplit(list.files(path="../raw/", pattern="_tot_cov.tab"), "_"), "[[", 1))
raw_pct <-  sort(do.call("c", lapply(all_chrs, getRawPct)))   

#### start loop here ###
for(ii in 1:30) 
	{
	#### data track blocks ####
	tdat       <- index[index$block %in% names(blocks[temp[ii,2][[1]]])] 
	nn         <- grep("^NN", colnames(mat)) 
	uc         <- grep("^UC", colnames(mat)) 
	tdat$nn    <- as.numeric(rowMeans(mat[names(tdat),nn]))  
	tdat$uc    <- as.numeric(rowMeans(mat[names(tdat),uc]))  
	tdat       <- tdat[,c("nn", "uc")]
	tbs        <- blocks[temp[ii,2][[1]]]
	all_exons  <- ucsc$full[ucsc$full$type=="exon"]

	#### data track full region  ####
	tpct       <- GRanges(paste0(unique(as.character(seqnames(tdat))) , ":", min(start(tdat)), "-", max(end(tdat))))
	tpct       <- raw_pct[queryHits(findOverlaps(raw_pct, tpct))]  
	nn         <- grep("^NN",  colnames(mcols(tpct))) 
	uc         <- grep("^UC",  colnames(mcols(tpct))) 
	tpct$nn    <- as.numeric(rowMeans(data.frame(mcols(tpct)[,nn]), na.rm=T))
	tpct$uc    <- as.numeric(rowMeans(data.frame(mcols(tpct)[,uc]), na.rm=T))
	tpct       <- tpct[, c("nn", "uc")]

	#### make sure the positions align with start end 
	raw_track        <- DataTrack(tpct, name="raw %") 
	meth_track       <- DataTrack(tdat, name="imp %")
	axis_track       <- GenomeAxisTrack()
	block_track      <- GeneRegionTrack(blocks[temp[ii,2][[1]]], name="blocks") 
	dmrseq_track     <- GeneRegionTrack(sig_dmrs[temp[ii,1]], name="dmrseq") 
	met_track        <- GeneRegionTrack(sig_mets[unique(queryHits(findOverlaps(sig_mets, tbs)))], name="metilene") 
	cpg_track        <- GeneRegionTrack(ucsc$cpg[unique(queryHits(findOverlaps(ucsc$cpg, tbs)))], name="cpg") 
	tfs_track        <- GeneRegionTrack(reg$TF_binding_site[unique(queryHits(findOverlaps(reg$TF_binding_site, tbs)))], name="tfs") 
	ctcf_track       <- GeneRegionTrack(reg$CTCF_binding_site[unique(queryHits(findOverlaps(reg$CTCF_binding_site, tbs)))], name="ctcf") 
	enh_track        <- GeneRegionTrack(reg$enhancer[unique(queryHits(findOverlaps(reg$enhancer, tbs)))], name="enhancers") 
	exon_track       <- GeneRegionTrack(all_exons[unique(queryHits(findOverlaps(all_exons, tbs)))], name="exon")


	displayPars(raw_track) <-  list("groups"=c("nn", "uc"), "col"=c("green", "red"), "col.frame"="lightgrey", "col.axis"="black") 
	displayPars(raw_track)  <- list("box.legend"=TRUE, "cex.legend"=1) 
	displayPars(meth_track) <- list("groups"=c("nn", "uc"), "col"=c("green", "red"), "col.frame"="lightgrey", "col.axis"="black") 
	displayPars(meth_track) <- list("legend"=FALSE) 

	all_tracks <- list(raw_track, meth_track, axis_track, block_track, dmrseq_track, met_track, cpg_track, tfs_track, ctcf_track, enh_track, exon_track)
	names(all_tracks) <- c("raw_track","meth_track","axis_track","block_track","dmrseq_track","met_track","cpg_track","tfs_track","ctcf_track","enh_track","exon_track")
	scheme     <- list("rotation.title"=0, "fontcolor.title"="black", "cex.title"=1,  "background.title" = "white", "frame"=TRUE) 

	for(i in 1:length(all_tracks)){displayPars(all_tracks[[i]]) <- scheme}
	displayPars(all_tracks$raw_track) <- list("rotation.title"=90) 
	displayPars(all_tracks$meth_track) <- list("rotation.title"=90) 
	displayPars(all_tracks$axis_track) <- list("frame"=FALSE) 
	displayPars(all_tracks$exon_track) <- list("stacking"="squish") 

	### remove empty tracks ###
	keep <- c(names(all_tracks)[1:3], names(which(lapply(all_tracks[4:length(all_tracks)], length) > 0)))
	all_tracks <- all_tracks[keep]

	pdf(paste0("../figures/test_", sprintf("%02d", ii), ".pdf"), width=10, height=10)  
	plotTracks(all_tracks, from=min(start(tdat))-5, to=max(end(tdat))+5, title.width=2)  
	dev.off() 
	#### end loop
	cat(ii, sep="\n") 
	}

### density plot ####
pdf("../figures/region_density.pdf", width=6, height=6) 
plot(density(width(blocks)), xlim=c(0,800), main="Size of regions", xlab="Width")  
lines(density(width(sig_dmrs)), col="red") 
lines(density(width(sig_mets)), col="blue") 
legend("topright", legend=c("blocks", "metiline", "dmrseq"), col=c("black", "blue", "red"), lty=1)
dev.off() 





