# swap fasta headers using two col dataframe 

#!/usr/bin/env Rscript

#changed shebang line to look for default Rscript loaded into environment


# this line is needed to run Rscript from command line
args <- commandArgs(trailingOnly = TRUE)

#	args[1] = shortname_fasta_file.faa
#	args[2] = keep_taxonomylookup.tab
#	args[3] = new longname_fasta_file.faa

		library(Biostrings)

		original_faa <- readDNAStringSet(args[1], format="fasta")
		short_to_long_query_genome_name_match <- read.table(args[2], header=FALSE, sep="\t")		
		names(original_faa) <- gsub("\\|.*","", names(original_faa))
	
		original_faa_headers <- as.data.frame(names(original_faa))
		#original_faa_headers <- data.frame(names(original_faa))
		colnames(original_faa_headers) <- c("V1")
		
		new_faa_headers <- merge(original_faa_headers, short_to_long_query_genome_name_match, by="V1", all.x=TRUE)
		new_faa_headers_sorted <- new_faa_headers[order(match(new_faa_headers[,1],original_faa_headers[,1])),]
		names(original_faa) <- new_faa_headers_sorted$V2
		writeXStringSet(original_faa, args[3], format="fasta", width=10000)
