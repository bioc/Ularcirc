% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ularcirc_standalone.R
\name{Junction_Sequence_from_Genome}
\alias{Junction_Sequence_from_Genome}
\title{Grab_BS_Junc_Sequence}
\usage{
Junction_Sequence_from_Genome(SelectUniqueJunct_Value, GeneList)
}
\arguments{
\item{GeneList}{: GeneList}

\item{SelectUniqueJunct_value}{: a dataframe with columns names startDonor, strandDonor, startAcceptor}
}
\description{
This function extracts genomic sequence that is likely to capture BSJ. Function does not cross
 validate to gene models.
}
