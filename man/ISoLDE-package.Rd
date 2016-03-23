\name{ISoLDE-package}

\alias{ISoLDE-package}
\alias{ISoLDE}

\docType{package}

\title{INTEGRATIVE STATISTICS OF ALLELE DEPENDENT EXPRESSION}

\description{
This package provides a new method for identifying genes with allelic bias.
This method is dedicated to data arising from RNA sequencing technologies.
The ISoLDE package implements original statistical methodology described in the
publication below.
}

\details{
ISoLDE method has been motivated by several literature limitations in taking
into account the data specificities and in making the most of biological
replicates. It is based on the definition of a new criterion using robust
estimation of the data variability. Variability estimation is of high
importance in statistical testing procedures because a difference significance
can only be assessed with regards to the intrinsic data variability. \cr

Two methods are available to identify allele specific expression: one is based
on bootstrap resampling while the second one uses an empirical threshold. The
first one is much more satisfying and is likely to give the most reliable
results but it can only be applied to data with at least three biological
replicates for each reciprocal cross. While strongly recommending to use at
least three replicates, the second method implements a robust solution when only
two replicates are available.
}

\author{
Christelle Reynès \email{christelle.reynes@igf.cnrs.fr},\cr
Marine Rohmer \email{marine.rohmer@mgx.cnrs.fr},\cr
Guilhem Kister \email{guilhem.kister@umontpellier.fr}
}

\references{
Reynès, C. et al. (2016): ISoLDE: a new method for identification of allelic
imbalance. \emph{Submitted}
}
