\name{isolde_test}
\alias{isolde_test}

\title{Statistical analysis of Allele specific read (ASR) counts}

\description{
  The main function of the ISoLDE package. Performs statistical test to identify
  genes with allelic bias and produces both graphical and textual outputs.
}

\usage{
  isolde_test(bias, method = "default", asr_counts, target,
             nboot = 5000, pcore = 75, graph = TRUE, ext = "pdf",
             text = TRUE, split_files = FALSE, prefix =
             "ISoLDE_result", outdir = "")
}

\arguments{
  \item{bias}{
    The kind of bias you want to study. It must be one of \dQuote{parental} or
    \dQuote{strain}.
  }
  \item{method}{
    specifies the statistical method to use for testing. It must be one of
    \dQuote{default} or \dQuote{threshold}. Default behaviour is to adapt to the
    number of replicates: when at least three biological replicates for each
    reciprocal cross are available the bootstrap resampling method is used, else
    the threshold method is applied. It is possible to force
    \code{\link{isolde_test}} to use the threshold method even when more than
    three replicates are available. In this case \code{method} must be set to
    \dQuote{threshold}. It is *not possible* to force a bootstrap method with
    less than three replicates. 
  }
  \item{asr_counts}{
    the \code{data.frame} containing the ASR counts to be tested.
    These data should be normalized and filtered (see the \code{\link{filterT}}
    function), although the function can run with non-normalized and
    non-filtered data (not recommended).
  }
  \item{target}{
    the target \code{data.frame} (obtained by the \code{\link{readTarget}}
    function).
  }
  \item{nboot}{
    specifies how many resampling steps to do for the bootstrap method.
    This option is not considered if \dQuote{threshold} value is set for
    \code{method}.
    Low values of \code{nboot} leads to less relevent results (default to 5000).
  }
  \item{pcore}{
    a value between 0 and 100 (default to 75) which specifies the proportion of
    cores (in percent) to be used for the bootstrap method.
   }
  \item{graph}{
    if \code{TRUE} (default) graphical outputs are produced (both on device and
    file). 
  }
  \item{ext}{
    specifies the extension of the graphical file output (does not work if
    graph = \code{FALSE}).
    It must be one of \dQuote{pdf} (default), \dQuote{png} or \dQuote{eps}.
  }
  \item{text}{
    if \code{TRUE} (default) textual output files are produced.
  }
  \item{split_files}{
    if text = \code{TRUE}, specifies if your want to have all genes in one same
    output file (\code{FALSE}, default) or four separate files according to the
    result: ASE, biallelic, undetermined or filtered (\code{TRUE}).
  }
  \item{prefix}{
    specifies the prefix for all output file names (default to "ISoLDE_result").
  }
  \item{outdir}{
    specifies the path where to write the output file(s) (default to current
    directory).
  }
}

\details{
  Before using this function, your data should be normalized and filtered
  (see the \code{\link{filterT}} function for filtering) although the function
  can run with non-normalized and/or non-filtered data.\cr

  The method depends on your minimum number of replicates for each reciprocal
  cross.\cr

  If only one replicate is found, the test can not be achieved and exits.\cr

  method=\dQuote{default} : If more than two replicates per cross, the method
  takes advantage of having enough   information by using bootstrap resampling
  to identify genes with allelic bias.\cr
  \cr
  If only two replicates are found in at least one cross, there is too few 
  information to obtain reliable distributions from resampling.
  Genes with allelic bias are identified thanks to empirically defined
  thresholds.\cr
  \cr
  method=\dQuote{threshold} : The empirical method will be processed instead of 
  the bootstrap one, even if more than two replicates per cross are found.\cr
  \cr
  Note that in differential RNA-seq analysis, at least three replicates are
  strongly recommended, as variability estimation quality is a key factor in
  statistical analysis.\cr
  \cr
  More details in Reynès, C. et al. (2016) ISoLDE: a new method for
  identification of allelic imbalance. \emph{Submitted}
}

\value{
  \item{listASE}{
    a \code{data.frame} with one row per gene (or transcript) identified as
    having an allelic bias and five columns: \cr
    - \dQuote{names} contains gene (or transcript) names such as
    \code{asr_counts} row\cr
    names, \cr
    - \dQuote{criterion} contains the criterion value (see vignette or Reynès
    et al. (2016)), \cr
    - \dQuote{diff_prop} the criterion numerator which contains the difference
    between proportions of either parents or strain origins, \cr
    - \dQuote{variability} the criterion denominator which quantifies the gene
    (or transcript) variability between replicates, \cr
    - \dQuote{origin} specifies the bias direction either "P" or "M" for
    parental bias or one of specified strain names for strain bias.
  }
  \item{listBA}{
    a \code{data.frame} with one row per gene (or transcript) identified as
    biallelically expressed and four columns corresponding to the first four
    ones in listASE. 	
  }
  \item{listUN}{
    a \code{data.frame} with one row per gene (or transcript) with undetermined
    status and six columns. The first five columns are the same as
    \code{listASE}, the last one may take three values: \cr
    - \dQuote{FLAG_consistency} for genes no statistical evidence of neither
    bias nor biallelic expression but whose parental or strain bias is always in
    the same direction across replicates, \cr
    - \dQuote{FLAG_significance} for genes with statistical evidence of bias but
    with discrepancies in bias direction across replicates, \cr
    - \dQuote{NO_FLAG} for other undetermined genes.
  }
  \item{listFILT}{
    a \code{data.frame} containing names of genes that have failed the minimal
    filtering step and thus that have not been considered during the statistical
    test.
  }
  ASE, BA and UN lists are sorted according to their criterion value. 
}

\references{
Reynès, C. et al. (2016): ISoLDE: a new method for identification of allelic
imbalance. \emph{Submitted}
}

\author{
  Christelle Reynès \email{christelle.reynes@igf.cnrs.fr},\cr
  Marine Rohmer \email{marine.rohmer@mgx.cnrs.fr}
}

\note{
  The bootstrap resampling step is performed many times (default to 5000).
  Hence, the function may run for a long time if performing the bootstrap method
  (until several minutes).
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
  data(filteredASRcounts)
  data(target)
  # Statistical analysis (forcing the threshold option)
  isolde_res <- isolde_test(bias = "parental", method = "threshold", 
asr_counts = filteredASRcounts, target = target, ext = "pdf",
prefix = "ISoLDE_test")
}
