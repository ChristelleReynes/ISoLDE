\name{filterT}
\alias{filterT}

\title{FILTERING GENES BEFORE STATISTICAL ANALYSIS}

\description{
Filter lowly expressed genes (or transcripts) according to a data driven
threshold, before any statistical
analysis. This step is not mandatory but strongly recommended.
}
\usage{
filterT(rawASRcounts, normASRcounts, target, tol_filter = 0,
bias)
}

\arguments{
\item{rawASRcounts}{
the \code{data.frame} containing raw counts (obtained with the
\code{\link{readRawInput}} function or any data.frame following
\code{\link{rawASRcounts}} format specifications).\cr
Raw count \code{data.frame} is required when filtering on raw or on normalized
data when the normalized data do not contain 0 counts.
(For simplicity purpose, we call '0 count' any value of zero in a count file).
}
\item{normASRcounts}{
the \code{data.frame} containing normalized counts (obtained with the\cr
\code{\link{readNormInput}} function or any data.frame following
\code{\link{normASRcounts}} format specifications).\cr
We strongly recommend to filter on normalized ASR counts.}

\item{target}{
the \code{data.frame} containing the target meta data (obtained with the\cr
\code{\link{readTarget}} function or any data.frame following
\code{\link{target}} format specifications).
}

\item{tol_filter}{
a value between 0 and 100 allowing to introduce tolerance rate into filtering
step:\cr
if tol_filter = 25 all genes having less than 25\% of their counts from at least
one parental (or strain) origin below the threshold are selected (the default
value 0 means all raw counts from at least one parental (or strain) origin must
be above threshold, 100 means that no filtering is applied).
}

\item{bias}{
The kind of allele expression bias you want to study. It must be one of
\dQuote{parental} or \dQuote{strain}.
}

}
\details{
Filtering in statistical analysis is recommended to avoid considering genes
(or transcript) without enough information, and thus to avoid a too strong
effect of multiple test correction.

The aim of our filtering method is to eliminate from analysis not enough
quantified genes, that is genes having mostly counts of 0 or near 0 for each
replicate in at least one condition (parent, strain). In this purpose, the
\code{\link{filterT}} function searches for the distribution of counts of a gene
in a condition when most of read counts are 0 for this condition.
This distribution allows to define a threshold. Hence, genes having less counts
than this threshold are eliminated.

The filtering step is not mandatory but strongly recommended. 
}

\value{
A list of two \code{data.frame}:
\item{\code{filteredASRcounts}}{This \code{data.frame} contains ASR counts that
have successfully passed the filtering step.}
\item{\code{removedASRcounts}}{This \code{data.frame} contains ASR counts that
have been removed by the filtering step.}\cr
\cr
Each line represents a feature (e.g. a gene, transcript).
Each column represents the number of allele-specific sens reads from either
the paternal or maternal parent for a given biological replicate, so that
you expect to have two columns per biological replicate.
}

\references{
Reynès, C. et al. (2016): ISoLDE: a new method for identification of allelic
imbalance. \emph{Submitted}
}

\author{
   Marine Rohmer \email{marine.rohmer@mgx.cnrs.fr},\cr
   Christelle Reynès \email{christelle.reynes@igf.cnrs.fr}
}

\note{
   \code{\link{filterT}} output on normalized data is the typical input for
\code{\link{isolde_test}}. 
}

\note{
A minimal filtering step will always be performed while applying the
\code{\link{isolde_test}} function.
It consists of eliminating all genes not satisfying these two conditions:\cr
- At least one of the two medians (of paternal or maternal ASR counts) is
different from 0;\cr
- There is at least one ASR count (different from 0) in each cross.
}

\examples{
# Loading all required data.frames
data(rawASRcounts)
data(normASRcounts)
data(target)

# Filtering genes from the ASR count data.frame in parental bias study
res_filterT <- filterT(rawASRcounts = rawASRcounts,
                       normASRcounts = normASRcounts,
                       target = target, bias="parental")
filteredASRcounts <- res_filterT$filteredASRcounts
removedASRcounts <- res_filterT$removedASRcounts
}


