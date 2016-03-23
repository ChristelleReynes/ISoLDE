\name{rawASRcounts}
\alias{rawASRcounts_file.txt}
\alias{rawASRcounts_file}
\alias{rawASRcounts}
\docType{data}

\title{RAW ALLELE SPECIFIC READ (ASR) COUNTS FOR AN EXPERIENCE WITH MORE THAN
TWO REPLICATES}

\description{
rawASRcounts_file.txt: A tab-delimited text file containing the raw values of
ASR counts for an experiment with more than two replicates.\cr
rawASRcounts.rda: the rawASRcounts_file.txt loaded into a data.frame by the
\code{\link{readRawInput}} function.
}

\format{
rawASRcounts_file.txt: A tab-delimited file.\cr
rawASRcounts.rda: A data.frame.
}

\details{
Each line represents a feature (e.g. a gene or a transcript).\cr
Each column represents the number of allele-specific sens reads from either\cr
the paternal or maternal parent for a given biological replicate, so that\cr
you expect to have two columns per biological replicate.\cr
Values in the matrix are raw allele-specific read counts.\cr
In case of double input, columns must be in the same order in both raw and
normalized ASR counts files.\cr
The rawASRcounts_file.txt file should be read and checked by the
\code{\link{readRawInput}} function.
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
\code{\link{readRawInput}} example of which uses the \code{\link{rawASRcounts}}
file.
}

\examples{
# rawASRcounts_file.txt
rawfile <- system.file("extdata", "rawASRcounts_file.txt", 
package = "ISoLDE")
rawASRcounts <- readRawInput(raw_file = rawfile, del = "tab", 
colnames = TRUE, rownames = TRUE)

# rawASRcounts.rda
data(rawASRcounts)
}

\keyword{datasets}
