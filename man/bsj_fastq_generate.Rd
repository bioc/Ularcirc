% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ularcirc_standalone.R
\name{bsj_fastq_generate}
\alias{bsj_fastq_generate}
\title{bsj_fastq_generate}
\usage{
bsj_fastq_generate(
  circRNA_Sequence,
  fragmentLength = 300,
  readLength = 100,
  variations = 4,
  headerID = ""
)
}
\arguments{
\item{circRNA_Sequence}{: Linear sequence of a circRNA. i.e. the backspice junction
is the first and last base of this sequence}

\item{fragmentLength}{: Is the length the library fragment}

\item{readLength}{: The sequence read length}

\item{variations}{: Number of sequences returned for each read type. Note each
sequence variation will start at a unique location (where possible)}

\item{headerID}{: Character identifier that will be incorporated into sequence header}
}
\value{
Returns a list of two DNAstring sets labelled "read1" and "read2" which correspond
to forward and reverse read pairs.
}
\description{
Takes a circRNA predicted sequence and generates synthetic short sequence reads
}
\examples{

library('Ularcirc')


# Generate a 500nt sequence containing A and which is flanked with GG and CC.
circRNA_Sequence <- paste(rep('A',500),collapse='')
circRNA_Sequence <- paste('GG',circRNA_Sequence, 'CC', sep='')
# The GG and CC ends of sequence represent ends of linear exons that are circularised.
# Therefore the backsplice junction (BSJ) is GGCC.
# Generate reads that alternate over this BSJ

fastqReads <- bsj_fastq_generate(circRNA_Sequence, fragmentLength=300, readLength=100,
               variations = 4,   # Four type I , II, III, and IV reads generated
               headerID='circRNA_example')  # Identifier incorporated in name of each sequence
# The following will indicate 12 sequences are present in each list entry
length(fastqReads$read1)
length(fastqReads$read2)

# Can create fastq file as follows
Biostrings::writeXStringSet( fastqReads$read1,"circRNA_Sample_R1.fastq.gz",
                              compress = TRUE, format="fastq")
Biostrings::writeXStringSet( fastqReads$read2,"circRNA_Sample_R2.fastq.gz",
                              compress = TRUE, format="fastq")
}
