## source("c:/Users/d.humphreys/OneDrive - Victor Chang/Documents/GitHub/Ularcirc_Bioconductor/Ularcirc/R/Ularcirc_miRNA.R")


######################################################################
#'
#' load_miRBase_fasta
#' This function extracts species specific mirbase sequences as a Biostrings RNAstringset. 
#'
#' @param mirbaselocalDir : Directory where miRbase files are located. If files are not detected will download requested files from mirBase
#' @param speciesCode : Three letter code of species (eg human = hsa). Used to subset species miRbase entries
#' @param filename : The miRBase filename. Expects this to be a fasta formatted file. If file does not exist will download from mirBase.
#'
#'
load_miRBase_fasta <- function(mirbaselocalDir = tempdir(), speciesCode = "hsa", filename="hairpin.fa" )
{
   # hairpinSeqs <- load_miRBase_fasta(mirbaselocalDir = tempdir(), speciesCode = "hsa", filename="hairpin.fa" ) 


	# Download hairpin.fa from miRBase
	url <- paste0("https://www.mirbase.org/download/", filename)
	local_file <- paste0(mirbaselocalDir, "/", filename)    # Local directory 
	if (! file.exists(local_file))   # does not exist
	{	test_connection <- tryCatch({ download.file(url, local_file, method = "auto") 
										TRUE
										},
										error = function(e) {
										FALSE
										})
										
		if (!test_connection) {
			warning(" Cannot reach miRBase server. Check your internet or if site is down try at a later time point.")
			return(NULL)
		}
	
	}	

	# Read the fasta file using Biostrings
	miRBase_sequences <- Biostrings::readRNAStringSet(local_file, format = "fasta")

	idx.species <- grep(paste0( "^", speciesCode), x=names(miRBase_sequences))
	
	# Subset for selected species
	return(miRBase_sequences[idx.species])
}




######################################################################
#'
#' circSequence
#' This function returns pre-defined full length circRNA sequences. 
#'
#' @param circRNA_ID : Character string of a circRNA_ID 
#' 
#'
#' @examples
#'       circSequence(circRNA_ID="slc8a1")
#'
#'
#' @export
#'
circSequence <- function(circRNA_ID="slc8a1")
{
    circRNAs <- list( slc8a1= "GAGAACATCTGGAGCTCGAGGAAATGTTATCGTTCCATATAAAACCATCGAAGGGACTGCCAGAGGTGGAGGGGAGGATTTTGAGGACACTTGTGGAGAGCTCGAATTCCAGAATGATGAAATTGT.TAGGTTGTGACAGTTGGAAGTGTCATGTACAACATGCGGCGATTAAGTCTTTCACCCACCTTTTCAATGGGATTTCATCTGTTAGTTACTGTGAGTCTCTTATTTTCCCATGTGGACCATGTAATT")
	
	idx <- which(names(circRNAs) %in% circRNA_ID)
	
	circRNASeq <- NULL

	if (length(idx) == 0)
	{   allIDs <- paste0(names(circRNAs), collapse=",")
		warning(paste0("circRNA ID not recognised. Try again using one of the following IDs: ", allIDs))
	}
	else
	{   circRNASeq <- circRNAs[idx]
	}

	return(circRNASeq)
}


######################################################################
#'
#' miR_binding_site_Analysis
#' This function analyses a genomic (RNA or DNA) sequence for matching miRNA seed sequences.
#' 
#'
#' @param Sequence_to_examine : RNA/DNA Sequence to be analysed for miRNA binding sites. Format can be a character string or Biostrings formatted object. 
#' @param speciesCode : Three letter code of species (eg human = hsa). Used to subset species miRbase entries
#' @param seed_length : Length of seed sequence of a miRNA (default 6).
#' @param seed_start : Starting position of seed sequence relative to mature miRNA sequence (default 2).
#' @param selected_miRs : names of miRNA.
#' @param rev_comp : Should function perform Reverse complement of Sequence_to_examine (default TRUE).
#' @param mirbaselocalDir : Directory where miRbase files are located. If files are not detected will download requested files from mirBase
#'
#' @details
#' 
#' The function will convert T to U and scan for matching miRNA seed sequences. 
#' Will return a list containing the following outputs: 
#'           SeedMatchResult  - counts of all seeds. Value of -1 indicates seed not detected. 
#'           Total_miR_number - the total number of miRNA scanned against. 
#'           miR_Seed_Lookup  - data frame of all miRNA seeds scanned against. 
#'
#' @examples
#'     Seq <- circSequence()
#'     output <- miR_binding_site_Analysis(Sequence_to_examine = Seq, species_code = "hsa", rev_comp = TRUE)
#'
#'     # display frequency of seed sequences.
#'     # Values of -1 do not occur
#'     head(output$SeedMatchResult)   
#'
#' @export
#'
miR_binding_site_Analysis <- function(Sequence_to_examine, species_code, seed_length=6, seed_start=2, selected_miRs = NULL, rev_comp = TRUE, mirbaselocalDir=tempdir())
{
	if (length(seed_length) == 0)
	{  seed_length <- 6 }
	
	stem_loops <- load_miRBase_fasta(mirbaselocalDir = mirbaselocalDir, speciesCode = species_code, filename="hairpin.fa" )
	mature_Seq <- load_miRBase_fasta(mirbaselocalDir = mirbaselocalDir, speciesCode = species_code, filename="mature.fa" )
	
	allseeds <- substr(x = mature_Seq, start = as.numeric(seed_start[1]), stop = as.numeric(seed_start[1])+ as.numeric( seed_length)-1)
	
	all_unique_seeds <- unique(allseeds)  #unique destroys names
    seed_df <- as.data.frame(allseeds)
	
	Sequence_to_examine <- Sequence_to_examine[[1]]  # Ensure there is only one entry
	Sequence_to_examine <- gsub(pattern = "T",replacement = "U",x = Sequence_to_examine)
	Sequence_to_examine <- Biostrings::RNAString(Sequence_to_examine)  
	
    if (rev_comp)
    { Sequence_to_examine <- reverseComplement(Sequence_to_examine)  }
    
    ## Find seed matches in sequence
    SeedMatchResult <-  mapply(FUN = function(x,y)
                    {   seed_regex <- paste("(?=",x,")",sep="")
                        return( gregexpr(seed_regex ,y, perl=TRUE) )
                    } , all_unique_seeds,as.character(Sequence_to_examine))

    return(list(SeedMatchResult=SeedMatchResult, Total_miR_number=length(all_unique_seeds), miR_Seed_Lookup=seed_df))
	
}
