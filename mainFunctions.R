
#-------------Sources----------#

source("elementeryFunctions.R")

#------------Functions---------#

baseCall_fun <- function(obj,fName,cb){
  
  Nucleotides<-c("A","C","G","T")
  abi_data <- NULL
  
  if (is.null(obj@data$DATA.9) ||
        is.null(obj@data$DATA.10) ||
        is.null(obj@data$DATA.11) ||
        is.null(obj@data$DATA.12) ||
        is.null(obj@data$PLOC.2) ||
        is.null(obj@data$PBAS.2) ||
        is.null(obj@data$PCON.2)){
    
    stop(paste(file,"is corrupted or doesn't contain RAW data !",sep=" "))
    
  }
  
  else{
    
    fName <- sub("([^.]+)\\.[[:alnum:]]+$", "\\1",fName)  
    
    tracematrix <- matrix(c(obj@data$DATA.9, 
                            obj@data$DATA.10, 
                            obj@data$DATA.11, 
                            obj@data$DATA.12), 
                          ncol=4)
    
    bcpos <- obj@data$PLOC.2 + 1
    
    #qual <- utf8ToInt(obj@data$PCON.2)
    qual <- obj@data$PCON.2
    
    rawSequence<-unlist(strsplit(obj@data$PBAS.2, ""))
    
    if (cb){
      Gpeaks <- getpeaks(tracematrix[,1])
      Apeaks <- getpeaks(tracematrix[,2])
      Tpeaks <- getpeaks(tracematrix[,3])
      Cpeaks <- getpeaks(tracematrix[,4])
      
      diffs <- diff(c(0,bcpos))
      
      starts <- bcpos - 0.5*diffs
      
      ends <- c(bcpos[1:(length(bcpos)-1)] + 
                  0.5*diffs[2:length(diffs)], 
                bcpos[length(diffs)] + 0.5*diffs[length(diffs)]
      )
      
      Sequence <- NULL
      
      for(i in 1:length(bcpos)) {
        
        if ( !rawSequence[i] %in% Nucleotides ){
          
          Apeak <- peakvalues(Apeaks, starts[i], ends[i])
          Cpeak <- peakvalues(Cpeaks, starts[i], ends[i])
          Gpeak <- peakvalues(Gpeaks, starts[i], ends[i])
          Tpeak <- peakvalues(Tpeaks, starts[i], ends[i])
          
          Aval <- maxvalues(tracematrix[,2], starts[i], ends[i])
          Cval <- maxvalues(tracematrix[,4], starts[i], ends[i])
          Gval <- maxvalues(tracematrix[,1], starts[i], ends[i])
          Tval <- maxvalues(tracematrix[,3], starts[i], ends[i])
          
          signals <- c(Apeak[1], Cpeak[1], Gpeak[1], Tpeak[1])
          maxes <- c(Aval[1], Cval[1], Gval[1], Tval[1])
          
          if(is.na(Apeak[2]) & is.na(Cpeak[2]) & is.na(Gpeak[2]) & is.na(Tpeak[2])){
            #Sequence<-c(Sequence," ")
            Sequence<-c(Sequence,baseReturn(which.max(maxes)))
            next 
          }
          
          
          if(max(signals)<100)
            #Sequence<-c(Sequence," ")
            Sequence<-c(Sequence,baseReturn(which.max(maxes)))
          else
            Sequence<-c(Sequence,baseReturn(which.max(signals)))
        }
        
        else
          
          Sequence<-c(Sequence,rawSequence[i])
      }
      sequenceOut <- Sequence
    }
    else
      sequenceOut <- rawSequence[-length(rawSequence)]
    
    abi_data <- list(fName, 
                     sequenceOut,
                     rawSequence[-length(rawSequence)],
                     tracematrix,
                     bcpos,
                     qual)
    names(abi_data) <- c("fileName", "seq", "seqraw", "tracematrix", "bcpos", "Q")
    
  }
  
  return (abi_data)
  
}

calcul.t5t3 <- function(qualVec , limit = 15 ){
  
  t3p<-t5p<-0
  
  for(i in 1:(length(qualVec) - 20)){
    vec <- qualVec[i : (i + 19)]
    if(length(which(vec < limit)) < 5){
      t5p <- i
      break
    }
  }
  
  for(i in length(qualVec): 20){
    vec <- qualVec[i : (i - 19)]
    if(length(which(vec < limit)) < 5){
      t3p <- length(qualVec)-i
      break
    } 
  }
  
  res <- list(t5=t5p,t3=t3p)
  return (res)
}

blaster<-function(seq,name,program){
  
  switch(program,
         "BLASTn"= {head <- "http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE=Nucleotides&PROGRAM=blastn&QUERY=%3E"},
         "BLASTx"= {head <- "http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE=Translations&PROGRAM=blastx&QUERY=%3E"},
         "TBLASTx"= {head <- "http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE=Translations&PROGRAM=tblastx&QUERY=%3E"},
         "MegaBLAST"= {head <- "http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE=MegaBlast&QUERY=%3E"})

  query<-paste(head,
               name,
               "%0a",
               char2string(seq),sep="")
  return(query)
}

pairwiseAligner<-function(seq1,seq2,type){
  
  Align <-pairwiseAlignment(seq1, seq2,type=type)
  capAlign <- capture.output(writePairwiseAlignments(Align,block.width=60))
  alignList <- list(program = gsub("# " ,"",capAlign[2]),
                header = c(gsub("# " ,"",capAlign[c(7:12)])," ",gsub("# " ,"",capAlign[c(14:18)])),
                align = capAlign[22:(length(capAlign)-2)])
  
  return(alignList)
}

chromHeight <- function(obj, trim5=0, trim3=0, width=100, pixelsperrow=150, showtrim=FALSE) {
  traces <- obj$tracematrix
  basecalls1 <- obj$seq
  basecalls2 <- obj$seqraw
  aveposition <- obj$bcpos
  basecalls1 <- basecalls1[1:length(aveposition)]
  basecalls2 <- basecalls2[1:length(aveposition)]
  if(showtrim == FALSE) {
    if(trim5+trim3 > length(basecalls1)) basecalls1 <- ""
    else basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
    if(trim5+trim3 > length(basecalls2)) basecalls2 <- ""
    else basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
    aveposition <- aveposition[(1 + trim5):(length(aveposition) - trim3)] 
  }
  indexes <- 1:length(basecalls1)
  trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all false if not trimmed
  if (!is.null(trim3)) {
    traces <- traces[1:(min(max(aveposition, na.rm=TRUE) + 10, nrow(traces))),]
  }
  if (!is.null(trim5)) {
    offset <- max(c(1, aveposition[1] - 10))
    traces <- traces[offset:nrow(traces),]
    aveposition <- aveposition - (offset-1)
  }
  
  valuesperbase <- nrow(traces)/length(basecalls1)
  tracewidth <- width*valuesperbase
  breaks <- seq(1,nrow(traces), by=tracewidth) 
  
  numplots <- length(breaks)
  return(numplots*pixelsperrow)
}

chromatogram.plot<-  function(obj, trim5, trim3, showOrigin = FALSE,
                              width, height=1, cex.mtext=1, cex.base=1, ylim=3,
                              filename=NULL, showtrim, showhets=TRUE) {
  
  
  
  originalpar <- par(no.readonly=TRUE)
  traces <- obj$tracematrix
  basecalls1 <- obj$seq
  basecalls2 <- obj$seqraw
  aveposition <- obj$bcpos
  showcalls <- ifelse(showOrigin, "both", "primary")
  basecalls1 <- basecalls1[1:length(aveposition)] 
  basecalls2 <- basecalls2[1:length(aveposition)]
  
  if(showtrim == FALSE) {
    if(trim5+trim3 > length(basecalls1)) basecalls1 <- ""
    else basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
    if(trim5+trim3 > length(basecalls2)) basecalls2 <- ""
    else basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
    aveposition <- aveposition[(1 + trim5):(length(aveposition) - trim3)] 
  }
  indexes <- 1:length(basecalls1)
  trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all 
  #false if not trimmed
  if (!is.null(trim3)) {
    traces <- traces[1:(min(max(aveposition, na.rm=TRUE) + 10, 
                            nrow(traces))), ]
  }
  if (!is.null(trim5)) {
    offset <- max(c(1, aveposition[1] - 10))
    traces <- traces[offset:nrow(traces),]
    aveposition <- aveposition - (offset-1)
  }
  maxsignal <- apply(traces, 1, max)
  ylims <- c(0, quantile(maxsignal, .75)+ylim*IQR(maxsignal))           
  p <- c(0, aveposition, nrow(traces))
  midp <- diff(p)/2
  starts <- aveposition - midp[1:(length(midp)-1)]
  starthets <- starts
  starthets[basecalls1 == basecalls2] <- NA
  ends <- aveposition + midp[2:(length(midp))]
  endhets <- ends
  endhets[basecalls1 == basecalls2] <- NA
  starttrims <- starts
  starttrims[!trimmed] <- NA
  endtrims <- ends
  endtrims[!trimmed] <- NA
  
  colortranslate <- c(A="green", C="blue", G="black", T="red")
  colorvector1 <- unname(colortranslate[basecalls1])
  colorvector1[is.na(colorvector1)] <- "purple"
  colorvector2 <- unname(colortranslate[basecalls2])
  colorvector2[is.na(colorvector2)] <- "purple"
  
  valuesperbase <- nrow(traces)/length(basecalls1)
  tracewidth <- width*valuesperbase
  breaks <- seq(1,nrow(traces), by=tracewidth)
  numplots <- length(breaks)
  if(!is.null(filename)) pdf(filename, width=8.5, height=height*numplots)
  par(mar=c(2,2,2,1), mfrow=c(numplots, 1))
  basecallwarning1 = 0
  basecallwarning2 = 0
  j = 1
  
  for(i in breaks) {
    range <- aveposition >= i & aveposition < (i+tracewidth)
    starthet <- starthets[range] - tracewidth*(j-1)
    starthet[starthet < 0] <- 0
    endhet <- endhets[range] - tracewidth*(j-1)
    endhet[endhet > tracewidth] <- tracewidth
    lab1 <- basecalls1[range]
    lab2 <- basecalls2[range]
    pos <- aveposition[range] - tracewidth*(j-1)
    colors1 <- colorvector1[range]
    colors2 <- colorvector2[range]
    starttrim <- starttrims[range] - tracewidth*(j-1)
    endtrim <- endtrims[range] - tracewidth*(j-1)
    plotrange <- i:min(i+tracewidth, nrow(traces))
    plot(traces[plotrange,1], type='n', ylim=ylims, ylab="", xaxt="n", 
         bty="n", xlab="", yaxt="n", , xlim=c(1,tracewidth))
    if (showhets==TRUE) {
      rect(starthet, 0, endhet, ylims[2], col='#D5E3F7', border='#D5E3F7')
    }
    if (showtrim==TRUE) {
      rect(starttrim, 0, endtrim, ylims[2], col='red', border='transparent', 
           density=15)
    }
    lines(traces[plotrange,1], col="black")
    lines(traces[plotrange,2], col="green")
    lines(traces[plotrange,3], col="red")
    lines(traces[plotrange,4], col="blue")
    mtext(as.character(which(range)[1]), side=2, line=0, cex=cex.mtext)
    
    for(k in 1:length(lab1)) {
      if (showcalls=="primary" | showcalls=="both") {
        if (is.na(basecalls1[1]) & basecallwarning1==0) {
          warning("Primary basecalls missing")
          basecallwarning1 = 1
        } 
        else if (length(lab1) > 0) {   
          axis(side=3, at=pos[k], labels=lab1[k], col.axis=colors1[k], 
               family="mono", cex=cex.base, line=ifelse(showcalls=="both", 0, 
                                                        -1), tick=FALSE)
        }
      }
      if (showcalls=="secondary" | showcalls=="both") {
        if (is.na(basecalls2[1]) & basecallwarning2 == 0) {
          warning("Secondary basecalls missing")
          basecallwarning2 = 1
        } 
        else if (length(lab2) > 0) { 
          axis(side=3, at=pos[k], labels=lab2[k], col.axis=colors2[k], 
               family="mono", cex=cex.base, line=-1, tick=FALSE)
        }
      }
    }
    j = j + 1
  }
  if(!is.null(filename)) {
    dev.off()
    
  }
  else par(originalpar)
}
                                  
cleanrefseq <- function(string) {
  string <- toupper(string)
  string <- gsub("[^ACGT]", "", string, perl=TRUE)
  return(string)
}

writeFasta<-function (sequences, names, file.out, open = "w", nbchar = 60, 
                      as.string = FALSE) 
{
  outfile <- file(description = file.out, open = open)
  write.oneseq <- function(sequence, name, nbchar, as.string) {
    writeLines(paste(">", name, sep = ""), outfile)
    if (as.string) 
      sequence <- s2c(sequence)
    l <- length(sequence)
    q <- floor(l/nbchar)
    r <- l - nbchar * q
    if (q > 0) {
      sapply(seq_len(q), function(x) writeLines(char2string(sequence[(nbchar * 
                                                                (x - 1) + 1):(nbchar * x)]), outfile))
    }
    if (r > 0) {
      writeLines(char2string(sequence[(nbchar * q + 1):l]), outfile)
    }
  }
  if (!is.list(sequences)) {
    write.oneseq(sequence = sequences, name = names, nbchar = nbchar, 
                 as.string = as.string)
  }
  else {
    n.seq <- length(sequences)
    sapply(seq_len(n.seq), function(x) write.oneseq(sequence = as.character(sequences[[x]]), 
                                                    name = names[x], nbchar = nbchar, as.string = as.string))
  }
  close(outfile)
}

detailSeq <- function(obj){
  
  if (!is.null(obj)){
  
  return (list(
    Format = paste("Format", obj@header@abif, sep=" : "),
    TUBE.1 = paste("Sample ID", obj@data$TUBE.1, sep=" : "),
    SMPL.1 = paste("Sample Name", obj@data$SMPL.1, sep=" : "),
    #Plate Label
    MODL.1 = paste("Instr. Model", obj@data$MODL.1, sep=" : "),
    MCHN.1 = paste("Instr. Name", obj@data$MCHN.1, sep=" : "),
    RUND.1 = paste("Run Start Date", paste(obj@data$RUND.1$month,obj@data$RUND.1$day,obj@data$RUND.1$year,sep="/"), sep=" : "),
    RUNT.1 = paste("Run Start Time", paste(obj@data$RUNT.1$hour,obj@data$RUNT.1$minute,obj@data$RUNT.1$second,obj@data$RUNT.1$hsecond,sep=":"), sep=" : "),
    RUND.2 = paste("Run Stop Date", paste(obj@data$RUND.2$month,obj@data$RUND.2$day,obj@data$RUND.2$year,sep="/"), sep="   :   "),
    RUNT.2 = paste("Run Stop Time", paste(obj@data$RUNT.2$hour,obj@data$RUNT.2$minute,obj@data$RUNT.2$second,obj@data$RUNT.2$hsecond,sep=":"), sep=" : "),
    LANE.1 = paste("Lane", obj@data$LANE.1, sep=" : "),
    NLNE.1 = paste("Number of Lanes", obj@data$NLNE.1, sep=" : "),
    #Spacing
    SN = paste("Signal Strs", paste( "A =",obj@data$`S/N%.1`[2],"C =",obj@data$`S/N%.1`[4],"G =",obj@data$`S/N%.1`[1] ,"T =",obj@data$`S/N%.1`[3],sep = " "),sep=" : "),
    PDMF.1 = paste("Mobility", obj@data$PDMF.1, sep=" : ")
    ))
  }
}

multiProcess <- function(dataFrameInput,recall,trim,trim5,trim3){
  
  processed<-list()
  
  unprocessed<-NULL
  
  for (i in 1 : length(dataFrameInput$name)){
    
    name <- gsub(" ", "_", dataFrameInput$name[i])
    
    obj<-read.abif(dataFrameInput$datapath[i])
    
    if (is.null(obj@data$DATA.9) ||
          is.null(obj@data$DATA.10) ||
          is.null(obj@data$DATA.11) ||
          is.null(obj@data$DATA.12) ||
          is.null(obj@data$PLOC.2) ||
          is.null(obj@data$PBAS.2) ||
          is.null(obj@data$PCON.2)){
      
      unprocessed<-c(unprocessed,dataFrameInput$name[i])
      next
      
    }
    
    else {
      
      data<-baseCall_fun(obj,name,recall)
      
      if (trim)
        t5t3 <- calcul.t5t3(data$Q)
      else
        t5t3 <- list(t5=trim5, t3=trim3)
      
      if (sum( t5t3$t5 , t5t3$t3 ) > length(data$seq))
        procSeq <- ""
      else
        procSeq <-  data$seq[(t5t3$t5+1):(length(data$seq)-t5t3$t3)]
      
      processed[[name]] <- procSeq
      
    }
    
  }
  
  result <- list (processed=processed,unprocessed=unprocessed)
  
  return(result)
  
}