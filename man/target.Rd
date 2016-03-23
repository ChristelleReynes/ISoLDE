\name{target}
\alias{target_file.txt}
\alias{target_file}
\alias{target}
\docType{data}

\title{METADATA ABOUT THE ASR COUNT DATA.}

\format{
target_file.txt: A tab-delimited text file.\cr
target.rda: A data.drame.
}

\description{
target_file.txt: A tab-delimited file describing your input data (raw and /
or normalized allele specific read (ASR) count file(s)).\cr
Each line of the target file corresponds to a column of the
\code{\link{rawASRcounts}} and / or \cr
\code{\link{normASRcounts}} \code{data.frames}.
Lines of target file MUST be in the same order as the columns in ASR count data.
Each line contains four values, separated by a character (e.g. a tabulation) :
\code{samples}, \code{parent}, \code{strain} and \code{replicate} (see the
Details section for more information).\cr
target.rda: The target_file.txt file loaded into a data.frame by the
\code{link{readTarget}} function.
}

\details{
Details of the three columns : \cr
\code{sample} : the biological sample name. A same sample name has to appear
twice in the target file : one line for the maternal allele and one line for the
paternal allele.\cr
\code{allele} : the parental origin of the ASR count. Two possible values:
\code{maternal} or \code{paternal}.\cr
\code{strain} : the strain origin of the ASR count. Exactly two different values
have to be provided in the whole file.\cr
\cr
The first line of the target file has to contain these column names in the same
order.\cr
These metadata are required for both \code{\link{filterT}} and
\code{\link{isolde_test}} functions.\cr
\cr 
Factice example: \cr
sample,parent,strain\cr
samp1,maternal,str1\cr
samp1,paternal,str2\cr
samp2,maternal,str1\cr
samp2,paternal,str2\cr
samp3,maternal,str1\cr
samp3,paternal,str2\cr
samp4,maternal,str1\cr
samp4,paternal,str2\cr
}


\author{
Marine Rohmer \email{marine.rohmer@mgx.cnrs.fr},\cr
Christelle Reynès \email{christelle.reynes@igf.cnrs.fr}
}

\seealso{
\code{readTarget} is a function to load into a \code{data.frame} and check the
input target file.
}
