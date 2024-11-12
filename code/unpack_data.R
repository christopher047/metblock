### download files ####

downloadReferenceFiles <- function() 
	{
	### download reference files ###

	### gencode for genes and transcripts ####
	gencode    = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.primary_assembly.annotation.gtf.gz"
	dest_gen   = "../reference/gencode.v36.primary_assembly.annotation.gtf.gz"
	if(!file.exists(dest_gen)){ download.file(gencode, destfile=dest_gen) }

	### ensembl regulatory features ####
	regulatory = "http://ftp.ensembl.org/pub/release-104/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz"
	dest_reg   = "../reference/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz"
	if(!file.exists(dest_reg)) { download.file(regulatory, destfile=dest_reg)}
	}

unpackRawData <- function() 
	{
	require(R.utils)	
	### unpack raw data as github limit is 25 MB ###
	input_files  <- list.files("../zipped", pattern="*.gz", full.name=T) 
	output_files <- gsub("../zipped", "../raw", input_files) 
	output_files <- gsub("\\.gz$", "", output_files) 
	for(ii in 1:length(input_files))
        	{
		if(file.exists(output_files[ii])){next}	
        	test_input   <- input_files[ii]
        	test_output  <- output_files[ii]
        	gunzip(test_input, destname = test_output, overwrite = FALSE, remove = FALSE)
        	}

	}

downloadReferenceFiles()
unpackRawData() 




