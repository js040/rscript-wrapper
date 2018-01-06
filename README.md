# rscript-wrapper

Skeleton-ish example of a wrapper for Rscripts. Purporse is to allow generation of Rscripts that can be run at the command line using flags. Because R does not currently allow the designation of flags when running Rscript from Shell (e.g example.r --input $FILEIN --output $FILEOUT

Example working command:
      ./fasta_headers_swap.py --shortnamefasta data/example-short-name.fna --keeptaxonomylookup data/name-key.list --newlongnamefastafile data/example-long-name.fna --pathtoscripts /Applications/ResearchSoftware/rscript-wrapper