\name{readTarget}
\alias{readTarget}

\title{READ THE TARGET FILE}

\description{Checks and loads into a \code{data.frame} your target input file.}

\usage{
readTarget(target_file, asr_counts, del = "\\t")
}

\arguments{
  \item{target_file}{
    A character-delimited text input file, containing metadata about ASR counts
    files (see \code{\link{target_file.txt}}).
}
  \item{asr_counts}{
    The \code{data.frame} containing values of ASR counts (obtained either by
    the \code{\link{readRawInput}} or the \code{\link{readNormInput}} function).
    It is used to perform checks on compatibility with the target file. 
}
  \item{del}{
    Specifies the delimiter for the target input file, usually a semi-colon ";",
    a coma "," or a tabulation "\\t". (default : "\\t").
    Note : None of your data values must contain this delimiter (be specially
    careful in gene names).
  }
}

\details{
See \code{\link{target_file.txt}} for more details about the \code{target_file}
format. 
}

\value{
a \code{data.frame} containing the target.
}

\author{
Marine Rohmer \email{marine.rohmer@mgx.cnrs.fr},\cr
Christelle Reynès \email{christelle.reynes@igf.cnrs.fr}
}

\seealso{
\code{\link{target_file.txt}}: the metadata file on which to run the
\code{\link{readTarget}} function.
}

\examples{
  # Target input file
  targetfile <- system.file("extdata", "target_file.txt", 
  package = "ISoLDE")
  # The data.frame containing ASR counts is also required
  data(rawASRcounts)
  # Load into a data.frame and check the target file
  target <- readTarget(target_file = targetfile, 
  asr_counts = rawASRcounts, del = "\\t")
}
