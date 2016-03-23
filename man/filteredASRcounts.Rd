\name{filteredASRcounts}
\alias{filteredASRcounts}
\alias{filteredASRcounts.rda}
\docType{data}

\title{
   NORMALIZED AND FILTERED ALLELE SPECIFIC READ (ASR) COUNTS FOR AN EXPERIENCE
WITH MORE THAN TWO REPLICATES
}

\description{
   A \code{data.frame} containing the normalized and filtered values of allele
   specific read counts,
   for an experiment with more than two replicates. This \code{data.frame} is
   obtained with the
   \code{\link{filterT}} function run on the \code{\link{normASRcounts}}
   \code{data.frame}.
}

\format{
   A \code{data.frame}.
}

\details{
   This \code{data.frame} is obtained with the \code{\link{filterT}} function.
   Each line represents a feature (e.g. a gene, transcript).
   Each column represents the number of allele specific sens reads from
   either the paternal or maternal parent for a given biological replicate,
   so that you expect to have two columns per biological replicate.
   Values in the matrix are filtered and normalized (RLE method) ASR counts.
}

\source{
   Extract from Bouschet, T. et al. (2016): In vitro corticogenesis from
   embryonic stem cells recapitulates the in vivo epigenetic control of
   imprinted gene expression. \emph{Submitted}  
   Subset of 6062 genes (after filtering).
}

\references{
   Bouschet, T. et al. (2016): In vitro corticogenesis from embryonic stem cells
   recapitulates the in vivo epigenetic control of imprinted gene expression.
   \emph{Submitted}
}

\seealso{
   \code{\link{filterT}}: a function to filter the ASR counts and produce the
   \code{filteredASRcounts}\cr
   \code{data.frame}.\cr
   \code{\link{isolde_test}}: example of which uses
   \code{\link{filteredASRcounts}}.
}

\examples{
   data(filteredASRcounts)
}

\keyword{datasets}
