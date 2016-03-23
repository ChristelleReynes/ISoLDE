#### ISoLDE Package

#  Copyright 2015 Christelle Reynes, Marine Rohmer
#  This file is part of ISoLDE.
#  ISoLDE is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#  ISoLDE is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License along with
#  this program.
#  If not, see <http://www.gnu.org/licenses/>.

## -------------------------------------------------------
## PUBLIC FUNCTION : READING AND CHECKING THE TARGET FILE
## -------------------------------------------------------

readTarget <- function(target_file, asr_counts, del = "\t")
{
  if (missing(target_file)) {
    stop("ERROR: Missing mandatory argument target_file.")
  }
  if (missing(asr_counts)) {
    stop("ERROR: Missing mandatory argument asr_counts.")
  }
  # Handle different ways to specify tabulation
  if ((del == "tab") || (del == "tabulation")) {
    del <- "\t"
  }

  message("Reading target file...")
  # Reading
  target <- read.table(target_file, header = TRUE, sep = del)
  
  # Checking the number of columns	OK
  if (ncol(target) != 3) {
    stop("ERROR : target file must contain 3 columns: sample, parent and strain.")
  } else
  {
	if (!all(colnames(target)==c("sample","parent","strain"))) stop("ERROR : target file must contain the following columns, in the same order: sample, parent and strain.")
  }

  # Checking the number of lines	OK
  if(ncol(asr_counts) != nrow(target)) {
    if (ncol(asr_counts) > nrow(target)) {
      stop("ERROR: missing lines in the target files. You need one header line plus one descriptive line per column of your data input file. See the vignette for more information.")
    } else if(ncol(asr_counts) < nrow(target)) {
      stop("ERROR: too much lines in the target files. You need one header line plus one descriptive line per column of your data input file. See the vignette for more information.")
    }
  }

  # Checking the target columns
  checkTarget(target, as.factor(target$parent))


  message("Done")
  return(target)
}
