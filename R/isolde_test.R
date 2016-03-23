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
## PUBLIC FUNCTION : STATISTICAL TEST TO FIND BIALLELIC OR IMPRINTED GENES
## -----------------------------------------------------------------------

isolde_test <- function(bias,
                        method = "default",
                        asr_counts,
                        target,
                        nboot = 5000, 
			pcore = 75, 
			graph = TRUE,
                        ext = "pdf",
                        text = TRUE,
                        split_files = FALSE,
                        prefix = "ISoLDE_result",
                        outdir = ""
                        ) {

  # Checking some parameters before continuing
  if (missing(asr_counts)) {
    stop("ERROR: Missing mandatory argument asr_counts.")
  }
  if (missing(target)) {
    stop("ERROR: Missing mandatory argument target.")
  }
  if (missing(bias)) {
    stop("ERROR: Missing mandatory argument bias.")
  } else {
    if ((bias != "parental") && (bias != "strain")) {
       stop("ERROR : Unexpected value for argument bias. Possible values : \"parental\" or \"strain\".")
    }
  }

  if ((method != "default") && (method != "threshold")) { 
    stop("The specified method is unrecognised. Possible values are \"default\" or \"threshold\".")
  }

  if (graph == TRUE) {
    if ((ext != "eps") && (ext != "pdf") && (ext != "png")) {
      stop(paste("Unknown extension for graphical output:",
                 ext,
                 "\nExpected: eps, png or pdf.",
                 sep = " "))
    }
  }
  if (pcore <= 0) {
	stop("ERROR: at least one core has to be used, please reconsider the pcore argument value (a percentage between 0 and 100).") 
  } else if (pcore < 10) {
	message(paste("Warning: pcore is a percentage between 0 and 100, be sure you want to use ",pcore,"% of your available cores.",sep="")) 
  } else if (pcore > 100) {
        stop(paste("ERROR: pcore is a percentage between 0 and 100. Found value: ",pcore,".",sep=""))
  }
    
  # If no outdir specified
  if (outdir == ""){
	outdir=getwd()
  }
  # Creates the outdir if it does not exist
  if (file.exists(outdir) == FALSE){
    dir.create(outdir)
  }

  # Replace unexpected spaces if any
  gsub(" ", "_", prefix)

  # ---------------------------------------------------------------------------------------------------------

  ### 1- MINIMAL FILTERING ###
  asr_counts=as.matrix(asr_counts)
  
  if (bias == "parental") {
      message("ISoLDE is searching for allele-specific gene expression due to parental biases...")
      ind1 <- getIndices(target, "maternal")
      ind2 <- getIndices(target, "paternal")
	  orig1 <- "M"
	  orig2 <- "P"
  } else if (bias == "strain") {
     message("ISoLDE is searching for allele-specific gene expression due to strain biases...")
     strain1 <- getStrainValues(target)$strain1
     strain2 <- getStrainValues(target)$strain2
     ind1 <- getIndices(target, strain1)
     ind2 <- getIndices(target, strain2)
	 orig1 <- strain1
	 orig2 <- strain2
  }
	minnbrep <- min(unlist(lapply(ind1, length)))
	
  if (minnbrep == 1) {
    stop("Sorry, it is not possible to apply this method without at least two biological replicates per cross")
  }
  
  ## filtering
  appFilter <- function(vec, iM, iP, seuil, tol) {
	    (sum(vec[unlist(iM)] >= seuil)>=length(unlist(iM))*(100-tol)/100) | (sum(vec[unlist(iP)] >= seuil)>=length(unlist(iP))*(100-tol)/100)
	  }
  repfilter <- apply(asr_counts, 1, appFilter, ind1, ind2, 2,0)	
  asr_counts_fil=asr_counts[repfilter,]
  
  
  # ---------------------------------------------------------------------------------------------------------

  ### 2- TESTING ###
  

  # PREDEFINED THRESHOLD METHOD (2 REPLICATES ONLY)
  if ((method == "default" && minnbrep == 2) || (method == "threshold")) {

    # If user chose the threshold method although each cross has more than 2 replicates
    if (method == "threshold" && minnbrep > 2) {
      message("The predefined threshold method is about to proceed, but we strongly recommand the bootstrap method because ISoLDE detected more than 2 replicates in each cross.")
      message("Testing...")
    } else {
      message("Given the small number of replicates, predefined thresholds will be applied. More replicates are recommended. See the vignette for more details.")
      message("Testing...")
    }

    ### Computing criterion and elements required for the graphical output
    crit9glob <- apply(asr_counts_fil, 1, stat9nonparf, ind1, ind2)
    Tdiff9glob <- unlist(lapply(crit9glob, function(listloc) listloc$diff)) 
    Tdenom9glob <- unlist(lapply(crit9glob, function(listloc) listloc$denom))
    Tcrit9glob <- unlist(lapply(crit9glob, function(listloc) listloc$crit))
    ## List creation (empirical thresholds)
    ASEthresh <- 1.7712
    BAthresh <- 0.7133
    resfinal <- rep(NA,nrow(asr_counts))
    names(resfinal) <- rownames(asr_counts)  
    listASE <- names(which(Tcrit9glob > ASEthresh))
    resfinal[listASE]="ASE"
    listBA <- names(which(Tcrit9glob < BAthresh))
    resfinal[listBA]="BA"
    listUN <- names(which(Tcrit9glob <= ASEthresh & Tcrit9glob >= BAthresh))
    resfinal[listUN]="UN"
    resfinal[rownames(asr_counts[repfilter==FALSE,])]="FILT"
    listFILT <- rownames(asr_counts[repfilter==FALSE,])
    ## Filtering according to consistency in print orientation
    # for ASE
	## filter for imprint consistency
	TsensOK=rep(TRUE,length(listASE))
	for (i in 1:length(listASE))
	{
	  Tsigloc=asr_counts[listASE[i],unlist(ind1)]-asr_counts[listASE[i],unlist(ind2)]
	  if (length(unique(sign(Tsigloc[Tsigloc!=0])))!=1) TsensOK[i]=FALSE
	}
	listASEnoncohe<-listASE[!TsensOK]
        resfinal[listASEnoncohe]="fakeASE"
	listASE<-listASE[TsensOK]
	listUN<-c(listUN,listASEnoncohe)
  
    # for undetermined
    repIND=names(which(Tcrit9glob <= ASEthresh & Tcrit9glob >= BAthresh))
    Tflag=rep(FALSE,nrow(asr_counts))
    names(Tflag)=rownames(asr_counts)	
	for (i in repIND)
	{
	  Tsigloc=asr_counts[i,unlist(ind1)]-asr_counts[i,unlist(ind2)]
	  if (length(unique(sign(Tsigloc)))==1) Tflag[i]=TRUE
	}
  }
  message("Done")
  # END PREDEFINED THRESHOLD (2 REPLICATES ONLY)

  # ---------------------------------------------------------------------------------------------------------

  # BOOTSTRAP METHOD (MORE THAN 2 REPLICATES)
  if (method == "default" && minnbrep > 2) {

    message("According to the satisfying number of replicates the full method can be applied. See the vignette for more details")
    message("Testing...")
    message("This step can last for a few minutes. Please be patient :)")
    message("Note: the full method uses a bootstrap step which means some results might change from one test to another.")

    if (nboot < 2000) {
      stop(paste("ISoLDE can not process the full method with less than 2000"," bootstrap steps because the result would be irrelevant.","\nDefault value: 5000.",sep = ""))
    }
    if (nboot > 10000) {
      message(paste("Note: An important bootstrap step number can last for a long time. You specified :",nboot,". Default value: 5000.",sep=""))
    }

	TNCore = pcore # % sued cores
	TNboot = nboot # Bootstrap runs
	TwithRepeat = 1
	Tthresh = 0.6
	enrres=matrix(NA,nrow(asr_counts_fil),3)
	enrpvalASE=matrix(NA,nrow(asr_counts_fil),3)
	for (runloc in 1:3)
		{
			TNRNF = t(asr_counts_fil)
			nIndM = sum(unlist(lapply(ind1,length)))
			Tdenom9glob=matrix(100001.1:(100000.1+nrow(asr_counts_fil)))
			Tdiff9glob=matrix(200001.1:(200000.1+nrow(asr_counts_fil)))
			Tcrit9glob=matrix(300001.1:(300000.1+nrow(asr_counts_fil)))
			nDDPH0 = choose(length(ind1[[1]]),2)+choose(length(ind1[[2]]),2)
			TdistdiffpropH0 = matrix(1.1, ncol = nrow(asr_counts_fil), nrow = nDDPH0)
			THistodistdiffpropH0 = matrix(0.0, ncol = 100)
			Tdistdiffcumul = matrix(0.0, ncol = 101)
			Tvecprob = matrix(0.0, ncol = 100)
			Tenrcrit9boot = matrix(0.0,nrow=TNboot,ncol=nrow(asr_counts_fil)) 
			TmoylocH0 = matrix(0.0,ncol=nrow(asr_counts_fil))
			Tpvalboot = matrix(0.0,ncol=nrow(asr_counts_fil))    

			########################################################
			### H1 test
			########################################################

			rep <- .Call("IsoldeP1", as.integer(TNCore),TNRNF, nrow(asr_counts_fil), ind1, ind2, Tcrit9glob, Tdiff9glob, Tdenom9glob,
				 TdistdiffpropH0, THistodistdiffpropH0, Tvecprob, Tdistdiffcumul, as.integer(TNboot), Tenrcrit9boot,
				 TmoylocH0, Tpvalboot, as.integer(TwithRepeat), PACKAGE="ISoLDE")
			
			##results extraction
			Tpvalbootfdr=p.adjust(t(Tpvalboot),method="BH")
			Tlistsign=rownames(asr_counts_fil[Tpvalbootfdr<=0.05,])
			TlistsignNum =which(Tpvalbootfdr<=0.05)
			enrpvalASE[,runloc]=Tpvalbootfdr		

			## filter for imprint consistency
			TsensOK=rep(TRUE,length(Tlistsign))
			for (i in 1:length(Tlistsign))
			{
			  Tsigloc=asr_counts[Tlistsign[i],unlist(ind1)]-asr_counts[Tlistsign[i],unlist(ind2)]
			  if (length(unique(sign(Tsigloc[Tsigloc!=0])))!=1) TsensOK[i]=FALSE
			}
			tabresH1=rep(NA,nrow(asr_counts_fil))
			tabresH1[Tpvalbootfdr<=0.05]="ASE"
			tabresH1[TlistsignNum[!TsensOK]]="fakeASE"

			########################################################
			### H0 test
			########################################################

			### replicates distribution for H1
			enrdiffH1<-NULL
			nbreadH1=asr_counts[Tlistsign,]

			TNRH1=t(nbreadH1)
			TenrdiffH1 = matrix(0.0, ncol = nrow(nbreadH1), nrow = length(ind1[[1]])*length(ind2[[2]]))
			Tenrcrit9bootH1bis=matrix(0.0,ncol=nrow(asr_counts_fil),nrow = TNboot)
			Tnbreadtot = matrix(0.0, ncol = nrow(asr_counts_fil), nrow = length(ind1[[1]])+length(ind1[[2]]))
			TmoylocH0X = matrix(0.0,ncol=nrow(asr_counts_fil))
			TpvalbootH1bis = matrix(0.0,ncol=nrow(asr_counts_fil))
					
			rep <- .Call("IsoldeP2", as.integer(TNCore),TNRNF, nrow(asr_counts_fil), ind1, ind2, Tcrit9glob, TNRH1, nrow(nbreadH1),
				TenrdiffH1, Tenrcrit9bootH1bis, Tnbreadtot, Tthresh, as.integer(TNboot), TmoylocH0X,
				TpvalbootH1bis, PACKAGE="ISoLDE")

			TpvalbootH1bis0.8fdr=p.adjust(t(TpvalbootH1bis),method="BH")
			
			## results extraction
			tabresH0=rep(NA,nrow(asr_counts_fil))
			tabresH0[TpvalbootH1bis0.8fdr<=0.05]="BA"
			tabresfin=rep(NA,nrow(asr_counts_fil))
			tabresfin[is.na(tabresH0) & !is.na(tabresH1)]=tabresH1[is.na(tabresH0) & !is.na(tabresH1)]
			tabresfin[!is.na(tabresH0) & is.na(tabresH1)]=tabresH0[!is.na(tabresH0) & is.na(tabresH1)]
			tabresfin[is.na(tabresH0) & is.na(tabresH1)]="UN"
			tabresfin[tabresH0=="BA" & tabresH1=="ASE"]="UN"
			tabresfin[tabresH0=="BA" & tabresH1=="fakeASE"]="UN"
			tabresfin[is.na(tabresH0) & tabresH1=="fakeASE"]="fakeASE"
			enrres[,runloc]=tabresfin
			
			## housekeeping
			rm(Tenrcrit9boot)
			rm(Tenrcrit9bootH1bis)
		}
	##############################	
	### final decision
	##############################
	repinco=apply(enrres,1,function(vec) length(unique(vec)))
	repinco=which(repinco>1)
	if (length(repinco)>0)
	    {
		for (i in repinco)
			{
				tabloc=table(enrres[i,])
				if (any(tabloc==2)) enrres[i,]=rep(names(tabloc)[which.max(tabloc)],3) else enrres[i,]=rep("UN",3)
			}
	    }
	resfinal=rep(NA,nrow(asr_counts))
	names(resfinal)=rownames(asr_counts)
	rownames(enrres)=rownames(asr_counts_fil)
	resfinal[rownames(enrres)]=enrres[,1]	
	resfinal[rownames(asr_counts[repfilter==FALSE,])]="FILT"
	### flag for consistent undetermined genes
	repIND=which(resfinal=="UN")
	Tflag=rep(FALSE,length(resfinal))
	for (i in repIND)
	{
	  Tsigloc=asr_counts[i,unlist(ind1)]-asr_counts[i,unlist(ind2)]
	  if (length(unique(sign(Tsigloc)))==1) Tflag[i]=TRUE
	}
	listASE=rownames(asr_counts)[resfinal=="ASE"]
    	listBA=rownames(asr_counts)[resfinal=="BA"]
    	listUN=c(rownames(asr_counts)[resfinal=="UN"],rownames(asr_counts)[resfinal=="fakeASE"])
        listFILT <- rownames(asr_counts[repfilter==FALSE,])
	
	names(Tcrit9glob)=rownames(asr_counts_fil)
	names(Tdiff9glob)=rownames(asr_counts_fil)
	names(Tdenom9glob)=rownames(asr_counts_fil)
	names(Tflag)=rownames(asr_counts)

    message("Done")    
  }
	
  
  # END BOOTSTRAP METHOD (MORE THAN 2 REPLICATES)

  # ---------------------------------------------------------------------------------------------------------
    ## Common information for both graphical and textual outputs

      current_date=format(Sys.time(), "%m-%d-%Y_%H-%M-%S")

      ## Adding differences, origins and tags for allele specific expressed genes (ASE)
      listASEtot <- cbind(Tcrit9glob[listASE],Tdiff9glob[listASE],Tdenom9glob[listASE])
      listASEtot <- as.data.frame(listASEtot)
      orig <- listASEtot[,2] > 0
      origin <- rep(orig2, length(orig))
      origin[orig] <- orig1
      listASEtot <- cbind(listASEtot, origin)
      listASEtot <- cbind(listASE, listASEtot)
      colnames(listASEtot)=c("name","criterion","diff_prop","variability","origin")
      rownames(listASEtot)=NULL
      ordloc=order(listASEtot[,2],decreasing=TRUE)
      listASEtot=listASEtot[ordloc,]   
    
      ## Adding differences, origins and tags for unknown genes (UN)
      listUNtot <- cbind(Tcrit9glob[listUN],Tdiff9glob[listUN],Tdenom9glob[listUN])
      listUNtot <- as.data.frame(listUNtot)
      orig <- listUNtot[,2] > 0
      origin <- rep(orig2, length(orig))
      origin[orig] <- orig1
      listUNtot <- cbind(listUNtot, origin)
      listUNtot <- cbind(listUN, listUNtot)
      listUNtot=cbind(listUNtot,factor(rep("NO_FLAG",nrow(listUNtot)),levels=c("NO_FLAG","FLAG_consistency","FLAG_significance")))
      listUNtot[names(which(Tflag==TRUE)),ncol(listUNtot)]="FLAG_consistency"
      listUNtot[rownames(asr_counts)[resfinal=="fakeASE"],ncol(listUNtot)]="FLAG_significance"
      rownames(listUNtot)=NULL
      colnames(listUNtot)=c("name","criterion","diff_prop","variability","origin","flag")
      ordloc=order(listUNtot[,2],decreasing=TRUE)
      listUNtot=listUNtot[ordloc,]


      ## Adding differences, origins and tags for unknown genes (BA)
      listBAtot <- cbind(Tcrit9glob[listBA], Tdiff9glob[listBA], Tdenom9glob[listBA])
      listBAtot <- as.data.frame(listBAtot)
      listBAtot <- cbind(listBA, listBAtot)
      colnames(listBAtot)=c("name","criterion","diff_prop","variability")
      rownames(listBAtot)=NULL
      ordloc <- order(listBAtot[,2], decreasing = TRUE)
      listBAtot <- listBAtot[ordloc,]
	  
      #Differenciate maternal and paternal ASE genes (for no difference : diff9_glob[listASE], denom9_glob[listASE])
      listASEtot_mat <- subset(listASEtot, listASEtot$origin == orig1)
      listASEmat <- listASEtot_mat[,1]
      listASEmat <- as.character(listASEmat)
      listASEtot_pat <- subset(listASEtot, listASEtot$origin == orig2)
      listASEpat <- listASEtot_pat[,1]
      listASEpat <- as.character(listASEpat)
      listUNcohe <- as.character(listUNtot[listUNtot$flag!="NO_FLAG",1])

      length_listFILT <- length(listFILT)

  # ---------------------------------------------------------------------------------------------------------

  
  
  ### 3- GRAPHICAL OUTPUT ###
  if (graph == TRUE) {
    plot(Tdiff9glob,
    Tdenom9glob,
    xlab = "allelic bias",
    ylab = "variability",
    type = "n")

    if ((minnbrep == 2) || (method == "threshold")) {
      seqdif <- seq(-1, 1, length = 500)
      seqdenom <- seq(0,0.5,length = 500)
      contour(seqdif,
              seqdenom,
              outer(seqdif, seqdenom^(1/2), "/"),
              levels = c(-ASEthresh, ASEthresh),
              col = "darkslategray",
              add = TRUE,lty=2)
      contour(seqdif,
              seqdenom,
              outer(seqdif, seqdenom^(1/2), "/"),
              levels = c(-BAthresh, BAthresh),
              col = "darkslategray",
              add = TRUE,lty=2)
    }    
    if (length(listASEmat)!=0) {
    points(Tdiff9glob[listASEmat],
           Tdenom9glob[listASEmat],
           col = "red3",
           pch = "+")
        if(length(listASEmat)<=60) {   
	   text(Tdiff9glob[listASEmat],
                Tdenom9glob[listASEmat],
                col = "red3",
		labels = listASEmat,
		pos = 1,
                cex = 0.6)
          }
    }
    if (length(listASEpat)!=0){   	   
       points(Tdiff9glob[listASEpat],
              Tdenom9glob[listASEpat],
	      col = "deepskyblue4",
              pch = "+") 
       if(length(listASEpat)<=60) {
          text(Tdiff9glob[listASEpat],
               Tdenom9glob[listASEpat],
               labels = listASEpat,
               pos = 1,
               col = "deepskyblue4",
               cex = 0.6)
       }
    }
    if (length(listBA)!=0){
       points(Tdiff9glob[listBA],
              Tdenom9glob[listBA],
              col = "darkmagenta",
              pch = "+")
    }
    if (length(listUN)!=0){
       points(Tdiff9glob[listUN[!listUN%in%listUNcohe]],
              Tdenom9glob[listUN[!listUN%in%listUNcohe]],
              col = "gray50",
              pch = "+")
    }
    if (length(listUNcohe)!=0){
       points(Tdiff9glob[listUNcohe],
              Tdenom9glob[listUNcohe],
              col = "gray50",
              pch = 10,cex=0.8)
    }	
    if (ext == "eps") {
       dev.copy2eps(file = paste(outdir,
                                 "/",
                                 prefix,
                                 "_graphical_output_",
                                 current_date,
                                 ".eps",
                                 sep = ""))
    } else if (ext == "png") {
       dev.copy(png, width=1200,
                height=1200,
                filename = paste(outdir,
                                 "/",
                                 prefix,
                                 "_graphical_output_",
                                 current_date,
                                 ".png",
                                 sep = ""))
    } else if (ext == "pdf") {
       dev.copy(pdf, file = paste(outdir,
                                  "/",
                                  prefix,
                                  "_graphical_output_",
                                  current_date,
                                  ".pdf",
                                  sep = ""))
    } else {
      stop(paste("Unknown extension for graphical output:",
           ext,
           "\nExpected: eps, png or pdf.",
           sep = " "))
    }
    dev.off()
  }

  #---------------------------------------------------------------------------------------------------------

  ### 4- TEXTUAL OUTPUT ###
  if (text == TRUE) {
     # Different files according to result
     if (split_files == TRUE) {   
	    # Allele specific expressed genes (ASE) textual output
	    if (!is.null(listASEtot)){
	       write.table(listASEtot,
		           file = paste(outdir, "/", prefix,"_ASE_", current_date, ".tsv", sep = ""),
		           quote = FALSE,
		           row.names = FALSE,
		           sep = "\t")
	    }

	    # Biallelic genes (BA) textual output
	    if (!is.null(listBAtot)){
	       write.table(listBAtot,
		          file = paste(outdir, "/", prefix, "_BA_", current_date, ".tsv", sep = ""),
		          quote = FALSE,
		          #col.names = "BA_genes",
		          row.names = FALSE,
		          sep = "\t")
	    }

	    # Undetermined genes (UN) textual output
	    if (!is.null(listUNtot)){
	       write.table(listUNtot,
		          file = paste(outdir, "/", prefix, "_UN_", current_date, ".tsv", sep = ""),
		          quote = FALSE,
		          row.names = FALSE,
		          sep = "\t")
	    }
	    
	    # Filtered genes (FILT) textual output
	    if (!is.null(listFILT)){
	       write.table(matrix(listFILT),
		          file = paste(outdir, "/", prefix, "_FILT_", current_date, ".tsv", sep = ""),
		          quote = FALSE,
		          row.names = FALSE,
			  col.names=FALSE,
		          sep = "\t")
            }

     # All genes in a same global file
     } else {
            global <- matrix(ncol=7)
            colnames(global) <- c("name","criterion","diff_prop","variability", "status","origin","flag")
            if (!is.null(listASEtot)){
               global <- merge(global,listASEtot,all=TRUE,sort=FALSE)
            }
            if (!is.null(listBAtot)){
               global <- merge(global,listBAtot,all=TRUE,sort=FALSE)
            }
            if (!is.null(listUNtot)){
               global <- merge(global, listUNtot, all=TRUE, sort=FALSE)
            }
            if (!is.null(listFILT)){
               listFILT <- as.data.frame(listFILT)
               colnames(listFILT) <- "name"
               global <- merge(global,listFILT,all=TRUE,sort=FALSE)
            } 
            # First line is only NA
            global <- global[2:dim(global)[1],]
            # Reorder as the inital count file
            global <- global[match(rownames(asr_counts), global$name),]
	    # fakeASE are UN genes tagged as FLAG_significance
            resfinal[resfinal=="fakeASE"] <- "UN"
            global$status <- resfinal
            global <- global[,c("name","criterion","diff_prop","variability", "status","origin","flag")]

       	    write.table(global,
		        file = paste(outdir, "/", prefix, "_ALL_", current_date, ".tsv", sep = ""),
		        quote = FALSE,
		        row.names = FALSE,
			col.names=TRUE,
		        sep = "\t")
       }
    }

    # ---------------------------------------------------------------------------------------------------------

    # Messages for user
    ASE_message <- paste(length(listASE),
                  " allele specific expressed (ASE) genes found and written in file ",
                  prefix,
                  "_ASE_",
                  current_date,
                  ".tsv.",
                  sep = "")
    BA_message <- paste(length(listBA),
                  " biallelic (BA) genes found and written in file ",
                  prefix,
                  "_BA_",
                  current_date,
                  ".tsv.",
                  sep = "")
    UN_message <- paste(length(listUN),
                  " Undetermined (UN) genes found and written in file ",
                  prefix,
                  "_UN_",
                  current_date,
                  ".tsv.",
                  sep = "")
    UNcohe_message <- paste("Among these undetermined (UN) genes, ",
                      length(listUNcohe),
                      " genes are coherent.\n",
                      "These genes correspond to the ones with the flag column set on \"FLAG_consistency\" of \"FLAG_significance\" in the undetermined genes file. ",
                      "See the vignette for more information on flags.",
                      sep = "")
    FILT_message <- paste(length_listFILT," genes have been filtered due to few reads.",sep="")

    message(paste(ASE_message, BA_message, UN_message, UNcohe_message, FILT_message, sep = "\n"))
  
  return(list("listASE"=listASEtot,"listBA"=listBAtot,"listUN"=listUNtot,"listFILT"=listFILT))
}
