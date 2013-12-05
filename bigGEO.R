#Read big GEO SOFT files
#Author: Tai-Hsien Ou Yang, Columbia University

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


#by Tai-Hsien Ou Yang

bigGEO <- function(fname,sampleSize=100) {
  gsmlist <- list()
  #gpllist <- list()
  GSMcount <- 0
  writeLines('Parsing....')
  con <- fileOpen(fname)
  lineCounts <- filegrep(con,"\\^(SAMPLE)",chunksize=10000)
  message(sprintf("Found %d entities...",nrow(lineCounts)))
  close(con)
  ## I close and reopen the file because on Windows, the seek
  ## function is pretty much guaranteed to NOT work
  con <- fileOpen(fname)
  ## This gets the header information for the GSE
  message(sprintf("Parsing Metadata..."))
  a <- readLines(con,lineCounts[1,1]-1) #-1 to avoid warning message
  header=parseGeoMeta(a)

  save(header,file="header.tmp")

  ## parse the actual entities, now
  for(j in 1:nrow(lineCounts)) {
    tmp <- strsplit(as.character(lineCounts[j,2])," = ")[[1]]
    accession <- tmp[2]
    message(sprintf("%s (%d of %d entities)",accession,j,nrow(lineCounts)))
    entityType <- tolower(sub("\\^","",tmp[1]))
    nLinesToRead <- lineCounts[j+1,1]-lineCounts[j,1]-1


    if(j==nrow(lineCounts))
      nLinesToRead <- NULL
    

      GSMcount=GSMcount+1
      #if(is.null(GSElimits)) {
      if( GSMcount %% sampleSize!=0){ 
        gsmlist[[accession]] <- .parseGSMWithLimits(con,n=nLinesToRead)
        
        }else{

        message(sprintf("Saving %d parsed samples...", GSMcount))
        save( gsmlist, file=paste( GSMcount,".tmp",sep="")  )
        rm(gsmlist)
        gc()
        gsmlist <- list()
        gsmlist[[accession]] <- .parseGSMWithLimits(con,n=nLinesToRead)
        }


   if(is.null(nLinesToRead)==TRUE){
      message(sprintf("Saving %d parsed samples...", GSMcount))
      save( gsmlist, file=paste( GSMcount,".tmp",sep="")  )
    }

  }
  close(con)

}





#Summarization by Tai-Hsien Ou Yang
summarizeGEO<-function(){

ptm <- proc.time()

data.sources = list.files(pattern="*.tmp")

load("header.tmp")
sample_id<-header$sample_id

#load( data.sources[1] ) #GPL platform data as the first
#geneSymbol<-make.unique(as.character(Table(gsmlist[[ 1 ]]@dataTable)[,11]))

load( data.sources[1] )

if( "ID" %in% names(Table(gsmlist[[1]]@dataTable))==TRUE ){
  probeset_ID<-as.character(Table(gsmlist[[1]]@dataTable)[,"ID"])
}else{
  probeset_ID<-as.character(Table(gsmlist[[1]]@dataTable)[,"ID_REF"])
}


#gsmlist[[]]@header$characteristics_ch1



  valName<- "VALUE"


ge<-matrix(NA,length(probeset_ID), 1  )
rownames(ge) <- probeset_ID



for( file.itr in 1:(length(data.sources)) ){ #1 is GPL platform data
  message(sprintf("Loading %d split...", file.itr ))

  load( data.sources[file.itr] )
  for( sample.itr in 1:length(gsmlist) ){
    if( valName  %in% names(Table(gsmlist[[ sample.itr ]]@dataTable))){ 
      ge.vec<-as.numeric(Table(gsmlist[[ sample.itr ]]@dataTable)[,valName])
      ge=cbind(ge, ge.vec)
      colnames(ge)[ncol(ge)]=names(gsmlist)[sample.itr]
    }
  }

}
ge=ge[,-1]


save(ge,file="ge.rda")

total.time =proc.time() - ptm
cat("takes",total.time[1],"sec\n")

}



annotateGEO<-function(ge){

  library("GEOquery")
  library("hgu133plus2.db")

  probeIDs<-ls(hgu133plus2SYMBOL)
  symbols<-mget(probeIDs,hgu133plus2SYMBOL)
  GPL570Symbols<-unlist(symbols)
  genes.shared=intersect(rownames(ge),names(GPL570Symbols))
  ge.reduced<-ge[genes.shared,]
  rownames(ge.reduced)<- GPL570Symbols[genes.shared]

  return(ge.reduced)
}




summarizeGEOClinical<-function(  survtext=c("follow_up_duration [(]years[)]: " ,"event_metastasis: "), survarrayid=c(5,6), clincaltext=NULL, clinicalarrayid=NULL, clinicalcol=NULL ){

ptm <- proc.time()

data.sources = list.files(pattern="*.tmp")

load("header.tmp")
sample_id<-header$sample_id

 # valName<- "VALUE"

surv<-matrix(NA,1, 2  )

if(is.null(clincaltext)==FALSE )
  clnc<-matrix(NA,1, length(clinicalarrayid) )


for( file.itr in 1:(length(data.sources)) ){ #1 is GPL platform data
  message(sprintf("Loading %d split...", file.itr ))

  load( data.sources[file.itr] )
  for( sample.itr in 1:length(gsmlist) ){
    #if( valName  %in% names(Table(gsmlist[[ sample.itr ]]@dataTable))){ 
      #ge.vec<-as.numeric(Table(gsmlist[[ sample.itr ]]@dataTable)[,valName])
      surv.1<-as.numeric(gsub(survtext[1], "", gsmlist[[sample.itr]]@header$characteristics_ch1[survarrayid[1]] ) )
      surv.2<-as.numeric(gsub(survtext[2], "", gsmlist[[sample.itr]]@header$characteristics_ch1[survarrayid[2]] ) )
               

      surv.vec<-matrix(c(surv.1,surv.2) ,1, 2  )
      surv=rbind(surv, surv.vec)

      #rownames(surv)[nrow(surv)]=strsplit(gsmlist[[sample.itr]]@header$title, "\\ ")[[1]][1]

      rownames(surv)[nrow(surv)]=names(gsmlist)[sample.itr]
      
      #if( is.na(surv.vec)[1] == TRUE)
      #  surv.vec<- c(as.numeric(gsub("drfs_even_time_years: ", "", gsmlist[[sample.itr]]@header$characteristics_ch1[14] ) ),
      #                              as.numeric(gsub("drfs_1_event_0_censored: ", "", gsmlist[[sample.itr]]@header$characteristics_ch1[13] ) ) )


      if(is.null(clincaltext)==FALSE ){

        #$characteristics_ch1
        #[1] "tissue: Breast Cancer"             "patient age: 51"                  
        #[3] "tumour size: 4"                    "nodes involved: 3"                
        #[5] "er status: 0"                      "tumour grade: 2"                  
        #[7] "distant-relapse event: 0"          "distant-relapse free survival: 10"
        clnc.vec=NULL
        for(i in 1:length(clinicalarrayid))
          clnc.vec<-c(clnc.vec,  as.numeric(gsub( clincaltext[i], "", gsmlist[[sample.itr]]@header$characteristics_ch1[clinicalarrayid[i]] ) ) )

          clnc=rbind(clnc, clnc.vec)
         # #rownames(clnc)[nrow(clnc)]=strsplit(gsmlist[[sample.itr]]@header$title, "\\ ")[[1]][1]

         rownames(clnc)[nrow(clnc)]=names(gsmlist)[sample.itr]
         if(is.null(clinicalcol)==F)
           colnames( clnc ) <- clinicalcol


      }
    
  }

}
surv=surv[-1,]
clnc=clnc[-1,]
#library("survival")
#surv<-Surv(surv)



save(surv,file="surv.rda")

if(is.null(clincaltext)==FALSE )
    save(clnc,file="clnc.rda")


total.time =proc.time() - ptm
cat("takes",total.time[1],"sec\n")

}


.example.summarizeGEOClinical<-function(){
summarizeGEOClinical( survtext=c( "relapse free survival time_days: ", 
                                  "relapse free survival event: "),
                      survarrayid=c(7,8),
                      clincaltext=c("er_status: ", "lymph node status: ", "size: " ),
                      clinicalarrayid=c(3, 6, 4),
                      clinicalcol=c("er","ln","size")
                        )
}



bigPipe<-function( fname, sampleSize=100 ){
  bigGEO(fname, sampleSize )
  summarizeGEO()
  load("ge.rda")
  ge<-annotateGEO(ge)
  return(ge)
} 

