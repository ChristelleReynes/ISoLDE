\name{readNormInput}
\alias{readNormInput}

\title{READ THE NORMALIZED DATA FILE}


\description{Checks and loads into a \code{data.frame} the input file containing
normalized allele-specific read (ASR) counts so that it can be input into
\code{\link{filterT}} and \code{\link{isolde_test}}.}


\usage{readNormInput(norm_file, del = "\t", rownames = TRUE, colnames = 
TRUE, dec = ".")}

\arguments{
  \item{norm_file}{
    A character-delimited input file containing normalized counts such as
    described in \code{\link{normASRcounts_file}}.
  }
  \item{del}{
    Specifies the delimiter for the input file, usually a semi-colon ";",
    a coma "," or a tabulation "\\t". (default : "\\t").
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
  \item{dec}{
    Specifies the character used in the file for decimal mark (default : ".").
  }
}

\value{
A \code{data.frame} containing normalized ASR counts from your input file.
}

\author{
Marine Rohmer \email{marine.rohmer@mgx.cnrs.fr},\cr
Christelle Reynès \email{christelle.reynes@igf.cnrs.fr}
}

\seealso{
\code{\link{normASRcounts_file.txt}}: the normalized ASR count file on which to
run the\cr
\code{\link{readNormInput}} function.\cr
\code{\link{readRawInput}}: a similar function for raw (non-normalized) ASR
count files.
}

\examples{
# character-delimited input file containing normalized ASR counts
normfile <- system.file("extdata", "normASRcounts_file.txt",
package = "ISoLDE")
# loading it into a data.frame using the readNormInput function
nbreadnorm <- readNormInput(norm_file = normfile, del = "tab", 
rownames = TRUE, colnames = TRUE, dec = ".") 
}
