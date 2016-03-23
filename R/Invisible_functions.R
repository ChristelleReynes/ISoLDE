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

## ********************************************************************************

## ---------------------------------------------------- ##
##                INVISIBLES FUNCTIONS                  ##
## ---------------------------------------------------- ##

## -------------------------------------------------
## INVISIBLE FUNCTION : TESTING ROW AND COLUMN NAMES

test_names <- function(cell, motif, value, problem, action, what) {
   if (length(grep(motif, cell, value=TRUE)) > 0) {
     if (action == "warning") {
       warning(paste("You specified ",what," = ",
                      value,
                      ", but ",
                      problem,
                      sep = ""))
     } else if (action == "stop") {
       stop(paste("You specified ",what," = ",
                      value,
                      ", but ",
                      problem,
                      sep = ""))
    }
  }
}

## ---------------------------------------
## INVISIBLE FUNCTION : GET CROSSES VALUES

getCrossValues <- function(target,cross) {
 #If we don't have 2 cross values
  if (length(levels(target$cross)) != 2) {
    values_found <- paste(levels(target$cross), collapse = " ")
    stop(paste("ERROR: cross column must contain 2 (and only 2) different values.",
               "Values found:",
               values_found,
               sep = "\n"))
  } else {
    return(list(cross1 = levels(target$cross)[1], cross2 = levels(target$cross)[2]))
  }
}

## ----------------------------------------
### INVISIBLE FUNCTION : GET STRAINS VALUES

getStrainValues <- function(target) {
  # If we don't have 2 strain values
  if (length(levels(target$strain)) != 2) {
    values_found <- paste(levels(target$strain), collapse = " ")
    stop(paste("ERROR: strain column must contain 2 (and only 2) different values.",
               "\nValues found:",
               values_found,
               sep = " "))
  } else {
    return(list(strain1 = levels(target$strain)[1], strain2 = levels(target$strain)[2]))
  }
}

## ----------------------------------------
## INVISIBLE FUNCTION : CHECK TARGET COLUMNS

checkTarget <- function(target, parent) {

  # generation of cross factor from the target file
  cross=NULL
  for (i in 1:nrow(target))
     {
	if (target$parent[i]=="maternal") {
           cross=c(cross,paste(target$strain[i],unique(target$strain[target$strain!=target$strain[i]]),sep="_"))
        } else { 
           cross=c(cross,paste(unique(target$strain[target$strain!=target$strain[i]]),target$strain[i],sep="_"))
        }
     }

  # Checking we have only two strains
  if (length(unique(target$strain))!=2) 
 	{
		stop(paste("ERROR: there should be exactly two strain names in the strain column."))
        }

  # Checking if we have the two strains for all maternal and all paternal
  if(length(unique(target$strain[target$parent=="maternal"])) < 2) {
	stop(paste("ERROR: Found the same strain name for all maternal lines which is not what is expected in an experiment involving reciprocal crosses. ",sep=""))
  }
  if(length(unique(target$strain[target$parent=="paternal"])) < 2) {
	stop(paste("ERROR: Found the same strain name for all paternal lines which is not what is expected in an experiment involving reciprocal crosses. ",sep=""))
  }

  
  # Combine the sample, parent and strain.
  replicate_names <- c()
  for (i in 1:nrow(target)) {
    replicate_names[length(replicate_names[]) + 1] <- paste(cross[i],
                                                            target[i,]$sample,
                                                            target[i,]$parent,
                                                            target[i,]$strain,
                                                            sep = ":")
  }
  replicate_names <- levels(as.factor(replicate_names))

# Checking that we have "maternal" and "paternal" for each replicate (and not another value)
  for (rep in replicate_names) {
    maternal <- FALSE
    paternal <- FALSE
    parent_tab <- c()
    for (i in 1:nrow(target)) {
      #if (paste(cross[i], target[i,]$replicate, sep = "-") == rep) {
      if (target[i,]$sample == unlist(strsplit(rep,":"))[2]) {
        # To be case insensitive
        current_parent <- tolower(target[i,]$parent)
        if (current_parent == "maternal") {
          maternal <- TRUE
        } else if (current_parent == "paternal") {
          paternal <- TRUE
        }
        parent_tab[length(parent_tab[]) + 1] <- current_parent
      }
    }
    values_found <- paste(levels(as.factor(parent_tab)), collapse = " ")
    if (!isTRUE(paternal)) {
      stop(paste("ERROR: missing parent \"paternal\" for sample : ",
                 unlist(strsplit(rep,":"))[2],
                 "\nFound: ",
                 values_found, sep=""))
    }
    if (!isTRUE(maternal)) {
      stop(paste("ERROR: missing parent \"maternal\" for sample : ",
                 unlist(strsplit(rep,":"))[2],
                 "\nFound: ", values_found, sep=""))
    }
    if (length(levels(as.factor(parent_tab))) != 2) {
      stop(paste("ERROR: For each sample, parent column of target file only requires both values: \"maternal\" and \"paternal\".\nDetails: ",
                 values_found, " values found for sample ",
                 unlist(strsplit(rep,":"))[1],
                 sep=""))
    }
    if (length(parent_tab) != 2) {
      stop(paste("ERROR: Parent column of target file requires two lines per sample, one for the \"maternal\" value and one for the \"paternal\" value.\nDetails: ",
                 length(parent_tab),
                 " values found for sample ",
                 unlist(strsplit(rep,":"))[2],
                 sep=""))
    }
  }
}


## --------------------------------
## INVISIBLE FUNCTION : GET INDICES
getIndices <- function(target, whichIndice)
{
  # Variables initialization
  strain1 <- getStrainValues(target)$strain1
  strain2 <- getStrainValues(target)$strain2
  
  if (whichIndice=="maternal")
     {
	indices=list(which(target$parent=="maternal" & target$strain==strain1),which(target$parent=="maternal" & target$strain==strain2))
     }
  if (whichIndice=="paternal")
     {
	indices=list(which(target$parent=="paternal" & target$strain==strain2),which(target$parent=="paternal" & target$strain==strain1))
     }
  if (whichIndice==strain1)
     {
	indices=list(which(target$parent=="maternal" & target$strain==strain1),which(target$parent=="paternal" & target$strain==strain1))
     }
  if (whichIndice==strain2)
     {
	indices=list(which(target$parent=="paternal" & target$strain==strain2),which(target$parent=="maternal" & target$strain==strain2))
     }
 return(indices)
}

## ----------------------------------------------
## INVISIBLE FUNCTION : CALCULATE TEST CRITERIONS

stat9nonparf <- function(vec, iM, iP) {
  diffp <- sum(vec[iM[[1]]]) / sum(vec[c(iM[[1]], iP[[1]])]) - sum(vec[iP[[2]]]) / sum(vec[c(iM[[2]], iP[[2]])])
  vecM <- vec[unlist(iM)]
  meM <- median(vecM)
  diffM <- abs(vecM - meM)
  if (sum(vecM) == 0) CVM <- 0 else CVM <- median(diffM) / sum(vecM)
  vecP <- vec[unlist(iP)]
  meP <- median(vecP)
  diffP <- abs(vecP - meP)
  if (sum(vecP) == 0) CVP <- 0 else CVP <- median(diffP) / sum(vecP)
  denom <- sum(c(CVM, CVP))
  return(list(diffp = diffp, denom = denom, crit = abs(diffp) / sqrt(denom)))
}

## -----------------------------------
## INVISIBLE FUNCTION : BOOTSTRAP STEP

repart <- function(val, probloc) {
  sum(sample(c(1,0), floor(val), replace = TRUE, prob = c(probloc, 1 - probloc)))
}

bootH0prob <- function(vec, distH0, vecbase) {
  probchoice <- sample(vecbase, 1, prob = distH0)
  aff1 <- apply(matrix(vec, ncol = 1), 1,repart, probloc = probchoice)
  aff2 <- floor(vec) - aff1 + (vec - floor(vec)) / 2
  aff1 <- aff1 + (vec - floor(vec))/2
  critloc <- stat9nonparf(c(aff1, aff2))$crit
  return(critloc)
}
