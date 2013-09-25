source("classes.R")

fileOpen <- function(fname,...) {
  con <- NULL
  if(length(grep('\\.gz$',fname,perl=TRUE))>0) {
    con <- gzfile(fname,open='rt')
  } else {
    con <- file(fname,'r')
  }
  return(con)
}


splitOnFirst <- function(x,pattern) {
  patlen <- nchar(pattern)
  matches <- regexpr(pattern,x)
  leftside <- substr(x,start=1,stop=matches-1)
  rightside <- substr(x,start=matches+patlen,stop=10000000)
  return(data.frame(leftside,rightside))
}

### This function does a grep on a file
### by doing a readline in chunks of size
### chunksize.
### Return value is a data.frame with the line number
### of each match and the line itself.
filegrep <-
  function(con,regex,chunksize=10000) {
    i <- 0
    ret <- NULL
    while(TRUE) {
      lines <- readLines(con,n=chunksize)


      if(length(lines)==0) {
        break
      }
      foundLines <- grep(regex,lines)
      foundTypes <- lines[foundLines]
      if(length(foundLines)>0) {
        foundLines <- foundLines+i
        tmp <- data.frame(foundLines=foundLines,foundTypes=foundTypes)
        if(is.null(ret)) {
          ret <- tmp
        } else {
          ret <- rbind(ret,tmp)
        }
      }
      i <- i+length(lines)
    }
    return(ret) 
  }

parseGeoColumns <- function(txt) {
  cols <- as.data.frame(splitOnFirst(txt[grep('^#',txt,perl=TRUE)],' = '))
  cols[,1] <- sub('#','',as.character(cols[,1]))
  colnames(cols) <- c('Column','Description')
  return(cols)
}

parseGeoMeta <- function(txt) {
  ## leader <- strsplit(grep('!\\w*_',txt,perl=TRUE,value=TRUE)[1],'_')[[1]][1]
  ## pull out only lines that are in the header
  tmp <- txt[grep("!\\w*?_",txt)]
  tmp <- gsub("!\\w*?_",'',tmp)
  first.eq <- regexpr(' = ',tmp)
  tmp <- cbind(substring(tmp,first=1,last=first.eq-1),
               substring(tmp,first=first.eq+3))
  tmp <- tmp[tmp[,1]!="",]
  header <- split(tmp[,2],tmp[,1])
  return(header)
}


fastTabRead <- function(con,sep="\t",header=TRUE,sampleRows=100,
                        colClasses=NULL,n=NULL,quote='"',...) {
### Need to read tables quickly, so guess the colclasses on the
### fly.  This is a bit dangerous since the first rows might not
### be representative of the entire set, but it is SO MUCH FASTER
### than the alternative, I have to do it.
  dat3 <- data.frame(NULL)
  numberOfLines <- -1
  if(!is.null(n)) {
    numberOfLines <- n-sampleRows
  }
  if(is.null(colClasses)) {
    if(!is.null(n)) {
      sampleRows <- min(sampleRows,n)
    }
    dat1 <- read.table(con,sep=sep,header=header,nrows=sampleRows,fill=TRUE,check.names=FALSE,
                       quote=quote,comment.char="",na.strings=c('NA','null','NULL','Null'),...)
    colclasses <- apply(dat1,2,class)
    colclasses[1] <- "factor"
    dat2 <- try(read.delim(con,sep=sep,colClasses=colclasses,
                       header=FALSE,quote=quote,comment.char="",
                       na.strings=c('NA','null','NULL','Null'),
                       nrows=numberOfLines,...),silent=TRUE)
    if(inherits(dat2,'try-error')) {
      dat3=dat1
    } else {
      colnames(dat2) <- colnames(dat1)
      dat3 <- rbind(dat1,dat2)
    }
  } else {
    dat3 <- read.delim(con,sep=sep,colClasses=colClasses,
                       header=header,quote=quote,comment.char="",
                       na.strings=c('NA','null','NULL',"Null"),nrows=numberOfLines,...)
  }
  return(dat3)
}


.parseGSMWithLimits <- function(con,n=NULL) {
  txt <- vector('character')
  i <- 0
  hasDataTable=FALSE
  while(i <- i+1) {
    tmp <- try(readLines(con,1))
    if(inherits(tmp,"try-error") | length(tmp)==0) {
      hasDataTable=FALSE
      break
    }
    txt[i] <- tmp
    if(length(grep('!\\w+_table_begin',txt[i],perl=TRUE))>0) {
      hasDataTable=TRUE
      break
    }
  }
  cols <- parseGeoColumns(txt)
  meta <- parseGeoMeta(txt)
  geoDataTable <- new("GEODataTable",columns=data.frame(),table=data.frame())
  if(hasDataTable) {
    nLinesToRead <- NULL
    if(!is.null(n)) {
      nLinesToRead <- n-length(txt)
    }
    dat3 <- fastTabRead(con,n=nLinesToRead)
    geoDataTable <- new('GEODataTable',columns=cols,table=dat3[1:(nrow(dat3)-1),])
  } 
  gsm <- new('GSM',
             header=meta,
             dataTable = geoDataTable)
}

bigGEO <- function(fname,sampleSize=100) {
  gsmlist <- list()
  #gpllist <- list()
  GSMcount <- 0
  writeLines('Parsing....')
  con <- fileOpen(fname)
  lineCounts <- filegrep(con,"\\^(SAMPLE|PLATFORM)",chunksize=10000)
  message(sprintf("Found %d entities...",nrow(lineCounts)))
  close(con)
  ## I close and reopen the file because on Windows, the seek
  ## function is pretty much guaranteed to NOT work
  con <- fileOpen(fname)
  ## This gets the header information for the GSE
  message(sprintf("Parsing Metadata..."))
  a <- readLines(con,lineCounts[1,1]-1)
  header=parseGeoMeta(a)

  save(header,file="header.rda")

  ## parse the actual entities, now
  for(j in 1:nrow(lineCounts)) {
    tmp <- strsplit(as.character(lineCounts[j,2])," = ")[[1]]
    accession <- tmp[2]
    message(sprintf("%s (%d of %d entities)",accession,j,nrow(lineCounts)))
    entityType <- tolower(sub("\\^","",tmp[1]))
    nLinesToRead <- lineCounts[j+1,1]-lineCounts[j,1]-1

    if(j==nrow(lineCounts)) {
      nLinesToRead <- NULL
      save( gsmlist, file=paste( GSMcount,".rda",sep="")  )
    }


    #if(entityType=="sample") {
    

      GSMcount=GSMcount+1
      #if(is.null(GSElimits)) {
      if( GSMcount %% sampleSize!=0){ 
        gsmlist[[accession]] <- .parseGSMWithLimits(con,n=nLinesToRead)
  
        }else{

        message(sprintf("Saving %d parsed samples...", GSMcount))
        save( gsmlist, file=paste( GSMcount,".rda",sep="")  )
        rm(gsmlist)
        gc()
        gsmlist <- list()
        gsmlist[[accession]] <- .parseGSMWithLimits(con,n=nLinesToRead)
        }



      #} else {

      #  message(sprintf("Not samples"))
      #  if((GSMcount>=GSElimits[1]) & (GSMcount<=GSElimits[2])) {
      #    gsmlist[[accession]] <- .parseGSMWithLimits(con,n=nLinesToRead)
      #  } else {
      #    if(!is.null(nLinesToRead)) {
      #      readLines(con,n=nLinesToRead)
      #    }
      #  }
      #}
    #}
    #if(entityType=="platform") {
      #gpllist[[accession]] <- .parseGPLWithLimits(con,n=nLinesToRead)
    #}
  }
  close(con)
  #return(new("GSE",
  #           header= header,
  #           gsms  = gsmlist,
  #           gpls  = gpllist))
}

