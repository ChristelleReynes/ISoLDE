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

## ---------------------------------------------------------
## PUBLIC FUNCTION : READING AND CHECKING THE RAW INPUT FILE
## ---------------------------------------------------------

readRawInput <- function(raw_file, del = "\t", rownames = TRUE, colnames = TRUE) {

  if (missing(raw_file)) {
    stop("ERROR: Missing mandatory argument raw_file.")
  }

  # Handle different ways to specify tabulation
  if ((del == "tab") || (del == "tabulation")) {
    del <- "\t"
  }
  rownames <- toupper(rownames)
  colnames <- toupper(colnames)

  # Reading
  message("Reading raw input file...")
  rawASRcounts <- ""
  if ((rownames == TRUE) && (colnames == TRUE)){
     rawASRcounts <- read.table(raw_file, sep = del, row.names = 1, header = TRUE)
     test_names(rownames(rawASRcounts)[1], "^[0-9]*$", rownames, "the first colum of your raw data file seems to contain counts (numeric) values.", "warning", "rownames")
     test_names(colnames(rawASRcounts)[1], "^X?[0-9]*$", colnames, "the first row of your raw data file seems to contain counts (numeric) values.", "warning", "colnames")
  } else if ((rownames == TRUE) && (colnames == FALSE)){
     rawASRcounts <- read.table(raw_file, sep = del, row.names = 1, header = FALSE)
     test_names(rawASRcounts[1,2], "[a-zA-Z]+", colnames, "your raw data file seems to also contain column names.", "stop", "colnames")
     test_names(rownames(rawASRcounts)[1], "^[0-9]$", rownames, "the first column of your raw data file seems to contain counts (numeric) values.", "warning", "rownames")
  } else if ((rownames == FALSE) && (colnames == FALSE)){
     rawASRcounts <- read.table(raw_file, sep = del, row.names = 0, header = FALSE)
     test_names(rawASRcounts[1,2], "[a-zA-Z]+", colnames, "your raw data file seems to contain column names.", "stop", "colnames")
     test_names(rawASRcounts[2,1], "[a-zA-Z]+", rownames, "your raw data file seems to contain row names.", "stop", "rownames")
  } else if ((rownames == FALSE) && (colnames == TRUE)){
     rawASRcounts <- read.table(raw_file, sep = del, row.names = 0, header = TRUE)
     test_names(rawASRcounts[2,1], "[a-zA-Z]+", rownames, "your raw data file seems to also contain row names.", "stop", "rownames")
     test_names(colnames(rawASRcounts)[1], "^X?[0-9]$", colnames, "the first row of your raw data file seems to contain counts (numeric) values.", "warning", "colnames")
  }

  # Checking the first lines of the data
  if (nrow(rawASRcounts) < 10) {
    lim = nrow(rawASRcounts)
  } else {
    lim = 10
  }
  for (i in 1:lim){
    for (j in 1:ncol(rawASRcounts)) {
      # Are values both digitals and integer ?
      if (length(grep("^[0-9]*$", rawASRcounts[i,j], value=TRUE)) == 0) {
        stop(paste("ERROR: Unexpected value in raw input file : \"",
                   rawASRcounts[i,j],
                   "\".\nOnly need integer values of raw ASRs counts (not normalized).",
                   "\nIf the first line or column of your data is supposed to be a header,",
                   "\nplease ensure it contains not only numbers.",
                   sep = ""))
      }
    }
  }
  # Conversion
  for (i in 1:ncol(rawASRcounts)) rawASRcounts[,i] <- as.numeric(as.character(rawASRcounts[,i]))
  rawASRcounts <- as.matrix(rawASRcounts)

  message("Done")
  return(rawASRcounts)
}
