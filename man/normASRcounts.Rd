\name{normASRcounts}
\alias{normASRcounts_file.txt}
\alias{normASRcounts_file}
\alias{normASRcounts}
\docType{data}

\title{NORMALIZED ALLELE SPECIFIC READ (ASR) COUNTS FOR AN EXPERIENCE WITH MORE
THAN TWO REPLICATES}

\description{
normASRcounts_file.txt: A tab-delimited text file containing the normalized
values of ASR counts for an experiment with more than two replicates.\cr
normASRcounts.rda: the normASRcounts_file.txt loaded into a data.frame by the\cr
\code{\link{readNormInput}} function.
}

\format{
normASRcounts_file.txt: A tab-delimited file.\cr
normASRcounts.rda: A data.frame.
}

\details{
Each line represents a feature (e.g. a gene or a transcript). \cr
Each column represents the number of allele-specific sens reads from either\cr
the paternal or maternal parent for a given biological replicate, so that\cr
you expect to have two columns per biological replicate.\cr
Values in the matrix are normalized (RLE method) ASR counts. \cr
In case of double input, columns must be in the same order in both raw and
normalized ASR counts files.\cr
The normASRcounts_file.txt file should be read and checked by the
\code{\link{readNormInput}} function.
}

\source{
Extract from Bouschet, T. et al. (2016): In vitro corticogenesis from embryonic
stem cells recapitulates the in vivo epigenetic control of imprinted gene
expression. \emph{Submitted} Subset of 6062 genes (after filtering).
}

\references{
Bouschet, T. et al. (2016): In vitro corticogenesis from embryonic stem cells
recapitulates the in vivo epigenetic control of imprinted gene expression.
\emph{Submitted}
}

\note{
A minimum of two biological replicates per cross is mandatory, however, we
strongly recommend to use more than two replicates per cross. This enables a
better estimation of variability and to use the bootstrap method to perform the
statistical test (see the \code{\link{isolde_test}} function).
}

\seealso{
\code{\link{readNormInput}} example of which uses the
\code{\link{normASRcounts}} file.
}

\examples{
# normASRcounts_file.txt
normfile <- system.file("extdata", "normASRcounts_file.txt",
package = "ISoLDE")
normASRcounts <- readNormInput(norm_file = normfile, del = "tab",
rownames = TRUE, colnames = TRUE)

# normASRcounts.rda
data(normASRcounts)
}

\keyword{datasets}
