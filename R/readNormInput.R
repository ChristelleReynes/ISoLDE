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

## ----------------------------------------------------------------
## PUBLIC FUNCTION : READING AND CHECKING THE NORMALIZED INPUT FILE
## ----------------------------------------------------------------

readNormInput <- function(norm_file, del = "\t", rownames = TRUE, colnames = TRUE, dec = ".") {

  if (missing(norm_file)) {
    stop("ERROR: Missing mandatory argument norm_file.")
  }
  # Handle different ways to specify parameters
  if ((del == "tab") || (del == "tabulation")) {
    del <- "\t"
  }
  rownames <- toupper(rownames)
  colnames <- toupper(colnames)

  # Reading 
  message("Reading normalized input file...")
  normASRcounts <- ""
  if ((rownames == TRUE) && (colnames == TRUE)){
     normASRcounts <- read.table(norm_file, sep = del, header = TRUE, row.names = 1, dec = dec)
     test_names(rownames(normASRcounts)[1], "^[0-9]*((\\.)[0-9]*)?$", rownames, "the first colum of your normalized data file seems to contain counts (numeric) values.", "warning", "rownames")
     test_names(colnames(normASRcounts)[1], "^X?[0-9]*((\\.)[0-9]*)?$", colnames, "the first row of your normalized data file seems to contain counts (numeric) values.", "warning"," colnames")
  } else if ((rownames == TRUE) && (colnames == FALSE)) {
     normASRcounts <- read.table(norm_file, sep = del, row.names = 1, dec = dec)
     test_names(normASRcounts[1,2], "[a-zA-Z]+", colnames, "your normalized data file seems to also contain column names.", "stop","colnames")
     test_names(rownames(normASRcounts)[1], "^[0-9]*((\\.)[0-9]*)?$", rownames, "the first column of your normalized data file seems to contain counts (numeric) values.", "warning", "rownames")
  } else if ((rownames == FALSE) && (colnames == TRUE)) {
    normASRcounts <- read.table(norm_file, sep = del, header = TRUE, dec = dec)
    test_names(normASRcounts[2,1], "[a-zA-Z]+", rownames, "your normalized data file seems to also contain row names.", "stop", "rownames")
    test_names(colnames(normASRcounts)[1], "^X?[0-9]*((\\.)[0-9]*)?$", colnames, "the first row of your normalized data file seems to contain counts (numeric) values.", "warning", "colnames")
  } else if ((rownames == FALSE) && (colnames == FALSE)) {
     normASRcounts <- read.table(norm_file, sep = del, header = FALSE, dec = dec)
     test_names(normASRcounts[1,2], "[a-zA-Z]+", colnames, "your normalized data file seems to contain column names.", "stop", "colnames")
     test_names(normASRcounts[2,1], "[a-zA-Z]+", rownames, "your normalized data file seems to contain row names.", "stop", "rownames")
  }

  # Checking the first lines
  if (nrow(normASRcounts) < 10) {
    lim <- nrow(normASRcounts)
  } else {
    lim <- 10
  }
  for (i in 1:lim) {
    for (j in 1:ncol(normASRcounts)) {
      # Error if coma found
      if (length(grep(".*,.*", normASRcounts[i,j], value = TRUE)) != 0) {
        stop(paste("ERROR: Unexpected coma in value:",
                   normASRcounts[i,j],
                   "\n Use the dot \".\" decimal character instead.",
                   sep = " "))
      }
      # Checking if there are some integer values
      # Some normalizations keep zeros as integer zeros
      if ((length(grep("^[0-9]*$", normASRcounts[i,j], value = TRUE)) != 0) && (grep("^[0-9]*$", normASRcounts[i,j], value = TRUE) != 0)) {
        warning(paste("Warning: Normalized input file seems to contain integer value:",
                       normASRcounts[i,j],
                       ".\nPlease ensure your data is normalized.",
                       sep = " "))
      }
      # Error if something else than numeric
      if (length(grep("^[0-9]*((\\.|,)[0-9]*)?$", normASRcounts[i,j], value = TRUE)) == 0) {
        stop(paste("ERROR: Found unexpected value in normalized input file:",
                   normASRcounts[i,j],
                   sep = " "))
      }
    }
  }
  # Conversion
  for (i in 1:ncol(normASRcounts)) {
    normASRcounts[,i] <- as.numeric(as.character(normASRcounts[,i]))
  }
  normASRcounts <- as.matrix(normASRcounts)
  normASRcounts <- normASRcounts 
  message("Done")
  return(normASRcounts)
}
