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

## -----------------------------------------------------------------------
## PUBLIC FUNCTION : DETERMINE THE THRESHOLD FOR FILTERING RAW ASR COUNTS
## -----------------------------------------------------------------------

filterT <- function(rawASRcounts=NULL, normASRcounts=NULL, target, tol_filter=0, bias)
{
  filter_on_raw = FALSE
  # Neither raw nor normalized data provided
  if (missing(rawASRcounts) & missing(normASRcounts)) {
    stop("ERROR: Missing mandatory argument rawASRcounts or normASRcounts.")
  # Filtering on normalized data without 0 counts (raw data needed)
  } else {
    if(length(grep("^[0-9]*((\\.)[0-9]*)$",rawASRcounts[1,1], value=TRUE)) > 0) {
       warning(paste("Warning: Data provided in rawASRcounts argument seems to contain normalized counts. Found non-integer value: ",rawASRcounts[1,1], sep=""))
    }
    if(length(grep("^[0-9]*$",normASRcounts[1,1], value=TRUE)) > 0) {
       warning(paste("Warning: Data provided in normASRcounts argument seems to contain raw counts. Found integer value: ",normASRcounts[1,1], sep=""))
    }
  }
  # Filtering on normalized data with 0 counts included (no raw data needed)
  if (missing(rawASRcounts) & !missing(normASRcounts)) {
    message("Note: Filtering on normalized data will only work if 0 counts are included.")
    rawASRcounts <- normASRcounts
    if(length(grep("^[0-9]*$", rawASRcounts[1,1], value=TRUE)) > 0) {
       warning(paste("Your data does not seem to be normalized. Found integer value: ", rawASRcounts[1,1], sep=""))
    }
  }
  # Filtering on raw data (no normalized data provided)
  if (missing(normASRcounts) & !missing(rawASRcounts)) {
    filter_on_raw = TRUE
    message("Note: Filtering on raw data. We strongly recommend to normalize your data before filtering it.")
  }
  if (missing(target)) {
    stop("ERROR: Missing mandatory argument target.")
  }
  if (!(is.numeric(tol_filter))) {
    stop(paste("ERROR: tol_filter must be a numeric value. Found: ", tol_filter, sep=""))
  }
  if ((tol_filter < 0) || (tol_filter > 100)) {
    print(tol_filter)
    stop("ERROR: tol_filter is a percentage between 0 and 100.")
  }
  if (tol_filter<1 & tol_filter!=0) {
	message(paste("Warning: tol_filter is a percentage between 0 and 100, be sure you want to tolerate ",tol_filter,"% of your data to be below threshold.",sep="")) 
  }
  if (tol_filter==100) {
	message(paste("Warning: tol_filter=100 corresponds to NO filtering step !",sep="")) 
  }
  if (missing(bias)) {
    stop("ERROR: Missing mandatory argument bias. Possible values: parental or strain.")
  } else if ((bias != "parental") && (bias != "strain")) {
    stop(paste("ERROR: bias argument must be one of parental or strain. Value found : ",bias, sep = ""))
  }


  message("Filtering your data...")
  raw_counts <- rawASRcounts

  # Get strain values
  strain1 <- getStrainValues(target)$strain1
  strain2 <- getStrainValues(target)$strain2

  # Get indices of maternals, paternals, strains1 and strains2
  indM <- getIndices(target, "maternal")
  indP <- getIndices(target, "paternal")
  indS1<- getIndices(target, strain1)
  indS2 <- getIndices(target, strain2)

  cpt0f<-function(vec) sum(vec == 0)
  cpt0M_tot <- apply(raw_counts, 1, cpt0f)
  reptout0 <- (cpt0M_tot != ncol(rawASRcounts))


  # max0 will contain the smallest count found in front of zeros
  max0 <- NULL
  raw_counts_fil <- raw_counts[reptout0,]
  appFilter <- function(vec, iM, iP, seuil, tol) {
	    (sum(vec[unlist(iM)] >= seuil)>=length(unlist(iM))*(100-tol)/100) | (sum(vec[unlist(iP)] >= seuil)>=length(unlist(iP))*(100-tol)/100)
	  }

  if (bias=="parental")
      {
	  for (i in 1:2) 
	     {
		# threshloc is a threshold depending on the replicate number.
		# We want each condition to have has at least 75% of 0 counts
		threshloc <- floor(0.75 * length(indM[[i]]))
		# searching for mothers
		cpt0M_1 <- apply(raw_counts_fil[, indM[[i]]], 1, cpt0f)
		# tmp_1 selects lines where there are enough 0 counts.
		tmp_1 <- raw_counts_fil[cpt0M_1 >= threshloc, indM[[i]]]
		max0 <- c(max0, apply(tmp_1, 1, max))
		# searching for fathers
		cpt0M_2 <- apply(raw_counts_fil[, indP[[i]]], 1, cpt0f)
		tmp_2 <- raw_counts_fil[cpt0M_2 >= threshloc, indP[[i]]]
		max0 <- c(max0, apply(tmp_2, 1, max))
  	     }
	  # Determining the final filtering threshold
	  threshold_filter=quantile(max0, 0.95)
          message(paste("Filtering threshold is ",threshold_filter,sep=""))
	  # Finally filtering normalized data using indices found on raw data
	  repfilter <- apply(rawASRcounts, 1, appFilter, indM, indP, threshold_filter, tol_filter)
      }
  else if (bias=="strain")
      {
	  for (i in 1:2) 
	     {
		# threshloc is a threshold depending on the replicate number.
		# We want each condition to have has at least 75% of 0 counts
		threshloc <- floor(0.75 * length(indS1[[i]]))
		# searching for strains 1
		cpt0M_1 <- apply(raw_counts_fil[, indS1[[i]]], 1, cpt0f)
		# tmp_1 selects lines where there are enough 0 counts.
		tmp_1 <- raw_counts_fil[cpt0M_1 >= threshloc, indS1[[i]]]
		max0 <- c(max0, apply(tmp_1, 1, max))
		# searching for strains 2
		cpt0M_2 <- apply(raw_counts_fil[, indS2[[i]]], 1, cpt0f)
		tmp_2 <- raw_counts_fil[cpt0M_2 >= threshloc, indS2[[i]]]
		max0 <- c(max0, apply(tmp_2, 1, max))
  	     }
	  # Determining the final filtering threshold
	  threshold_filter=quantile(max0, 0.95)
          message(paste("Filtering threshold is ",threshold_filter,sep=""))
	  # Finally filtering normalized data using indices found on raw data
	  repfilter <- apply(rawASRcounts, 1, appFilter, indS1, indS2, threshold_filter, tol_filter)	
      }
 
  if (filter_on_raw == FALSE) {
     filteredASRcounts <- normASRcounts[repfilter == TRUE,]
     removedASRcounts <- normASRcounts[repfilter == FALSE,]
  } else {
     filteredASRcounts <- rawASRcounts[repfilter == TRUE,]
     removedASRcounts <- rawASRcounts[repfilter == FALSE,]
  }

  message("Done")
  return(list("filteredASRcounts"=filteredASRcounts,"removedASRcounts"=removedASRcounts))
}
