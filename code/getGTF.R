
### get GTF ###
getGTFAnnotation <- function(gtf_file = gtf_file)
        {
        #this is hg38 version 36
        require(rtracklayer)
        require(GenomicRanges)
	gtf       <- readGFF(gzfile(gtf_file))
        gtf       <- GRanges(gtf)
        common    <- paste0("chr", c(1:22,"M", "X", "Y"))
        gtf       <- keepSeqlevels(gtf, common, pruning.mode="coarse")
        return(gtf)
        }

getGTF <- function(gtf_file="../reference/gencode.v36.primary_assembly.annotation.gtf.gz") 
	{
	require(GenomicRanges)
	require(TxDb.Hsapiens.UCSC.hg38.knownGene)
	require(org.Hs.eg.db)
	require(clusterProfiler) 
	gtf        <- getGTFAnnotation(gtf_file)

	### add entrez to simplify script ####
	all_genes  <- unique(gtf$gene_name)
	temp       <- bitr(all_genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)
	temp       <- setNames(temp$ENTREZID, temp$SYMBOL) 
	gtf$entrez <- "" 
        gtf$entrez <- as.character(temp[gtf$gene_name]) 

	### transcripts and genes ####
	full       <- gtf
	gxs        <- gtf[gtf$type == "gene"]
	gtf        <- gtf[gtf$type == "transcript"]
	names(gtf) <- gtf$transcript_id

	### load txdb ####
	txdb      <- TxDb.Hsapiens.UCSC.hg38.knownGene
	common    <- paste0("chr", c(1:22,"M", "X", "Y"))
	txdb      <- keepSeqlevels(txdb, common, pruning.mode="coarse")

	### sync txdb and trans ###
	trans        <- transcripts(txdb)
	names(trans) <- trans$tx_name
	gtf          <- gtf[intersect(names(gtf), names(trans))]
	trans        <- trans[intersect(names(trans), names(gtf))]
	trans        <- trans[names(gtf)]

	### check ####
	all(gxs$gene_name %in% gxs$gene_name)
	all(gtf$start == trans$start) 
	all(gtf$end   == trans$end) 

	### txdb ###
	prom2   <- promoters(gtf, upstream=2000, downstream=200)
	prom5   <- promoters(gtf, upstream=5000, downstream=200)
	prom10  <- promoters(gtf, upstream=10000, downstream=200)
	prom20  <- promoters(gtf, upstream=20000, downstream=200)
	
	### misc data ###
	tss     <- import.bed("../reference/refTSS_v3.1_human_coordinate.hg38.bed")
	
	### cpg #####
	cpg <- read.table("../reference/cpgIslandExtV36.txt", sep="\t")
	colnames(cpg) <- c("bin", "chr", "start", "end", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")
	cpg       <- GRanges(cpg)
	common    <- paste0("chr", c(1:22,"X", "Y"))
	cpg       <- keepSeqlevels(cpg, common, pruning.mode="coarse")
	
	### return ###
	x <- list("full"=full, "gtf"=gtf, "prom2"=prom2, "prom5"=prom5, "prom10"=prom10, "prom20"=prom20, "tss"=tss, "cpg"=cpg) 
	### return ###
	return(x)  
	}

getRegSites <- function(filename="../reference/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz") 
	{
	require(rtracklayer)
	reg <- import.gff(gzfile(filename))
	reg      <- keepSeqlevels(reg, c(1:22, "X", "Y"), pruning.mode="coarse")
	seqlevels(reg) <- paste0("chr", seqlevels(reg))
	reg <- sort(reg) 
	return(split(reg, as.factor(reg$type))) 
	}

