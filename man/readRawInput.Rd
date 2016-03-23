\name{readRawInput}
\alias{readRawInput}

\title{READ THE RAW DATA FILE}

\description{Checks and loads into a \code{data.frame} the input file containing
raw allele specific read (ASR) counts so that it can be input into
\code{\link{filterT}}.}

\usage{readRawInput(raw_file, del = "\t", rownames = TRUE, colnames =
    TRUE)}

\arguments{
  \item{raw_file}{
    A character-delimited input file containing raw ASR counts such as described
    in \code{\link{rawASRcounts_file}}.
  }
  \item{del}{
    Specifies the delimiter for the input file, usually a semi-colon ";", a coma
    "," or a tabulation "\\t". (default : "\\t").
    Note : None of your data values must contain this delimiter (be specially
    careful in gene names).
  }
  \item{rownames}{
    Specifies if the file contains some row names to consider.
    Possible values: TRUE or FALSE (default: TRUE).
  }
  \item{colnames}{
    Specifies if the file contains some column names to consider.
    Possible values: TRUE or FALSE (default: TRUE).
  }
}

\details{
Raw ASR counts are only required for the filtering step (with the
\code{\link{filterT}} function) in case the normalized data do not contain 0
counts anymore.
If you do not want to perform the filtering step or if you still have 0 counts
in your normalized file, you do not need to load raw ASR counts.
(For simplicity purpose, we call '0 count' any value of zero in a count file).
}

\value{
A \code{data.frame} containing raw ASR counts from your input file.
}

\author{
Marine Rohmer \email{marine.rohmer@mgx.cnrs.fr},\cr
Christelle Reynès \email{christelle.reynes@igf.cnrs.fr}
}

\seealso{
\code{\link{rawASRcounts_file.txt}}: the raw ASR count file on which to run the
\code{\link{readRawInput}} function. 
\code{\link{readNormInput}}: a similar function for normalized ASR count file.
}

\examples{
# character-delimited input file containing raw ASR counts
rawfile <- system.file("extdata", "rawASRcounts_file.txt", 
package = "ISoLDE")
# loading it into a data.frame using the readRawInput function
nbread <- readRawInput(raw_file = rawfile, del = "tab", 
rownames = TRUE, colnames = TRUE) 
}


