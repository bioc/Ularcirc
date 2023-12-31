% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/file_io.R
\name{chimericStats}
\alias{chimericStats}
\title{chimericStats}
\usage{
chimericStats(chimericDT)
}
\arguments{
\item{chimericDT}{: Data table of chimeric junctions as provided by
STAR aligner}
}
\description{
Simple function that returns a list of basic stats obtained from a STAR
chimeric file
}
\examples{
extdata_path <- system.file("extdata",package = "Ularcirc")
chimeric.file <- paste0(extdata_path,"/SRR444655_subset.Chimeric.out.junction.gz")
chimericDT <- Ularcirc::loadSTAR_chimeric(chimeric.file,returnColIdx = 1:14)
Ularcirc::chimericStats(chimericDT$data_set)
chimericDT$filtered <- Ularcirc::FilterChimericJuncs(chimericDT$data_set, canonicalJuncs = TRUE)
Ularcirc::chimericStats(chimericDT$filtered)


}
\seealso{
FilterChimericJuncs
}
