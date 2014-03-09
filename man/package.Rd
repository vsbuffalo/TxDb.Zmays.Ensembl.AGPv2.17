\name{TxDb.Zmays.Ensembl.AGPv2.17}
\docType{package}

\alias{TxDb.Zmays.Ensembl.AGPv2.17-package}
\alias{TxDb.Zmays.Ensembl.AGPv2.17}
\alias{TxDb.Zmays.Ensembl.AGPv2.17}


\title{Annotation package for TranscriptDb object(s)}

\description{
  This package loads one or more TranscriptDb objects.  Such TranscriptDb
  objects are an R interface to prefabricated databases contained by
  this package.

  The names of any objects exposed by this package indicate the origin and
  resources exposed.  So for example TxDb.Hsapiens.UCSC.hg19.knownGene
  would be a TranscriptDb object, of Homo Sapiens data from UCSC build
  hg19 based on the knownGene Track.
}

\author{Vince Buffalo}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.Zmays.Ensembl.AGPv2.17)
## list the contents that are loaded into memory
ls('package:TxDb.Zmays.Ensembl.AGPv2.17')
## show the db object that is loaded by calling it's name
TxDb.Zmays.Ensembl.AGPv2.17

}

\keyword{package}
\keyword{data}
