library(GenomicRanges)
library(regioneR)
library(dplyr)
maskCentromers <- function(GR, GR.centromers, extend.start = 0,
                           extend.end = 0){
  require(regioneR)
  if (extend.start != 0 & extend.end != 0){
    GR.centromers <- extendRegions(GR.centromers,
                                   extend.start=extend.start,
                                   extend.end = extend.end)
  }
  GR <- GR[!GR %over% GR.centromers]
  return(GR)
}
.chrAsNum <- function(tbl){
  tbl$chrom <- gsub("chr", "", tbl$chrom)
  tbl$chrom[tbl$chrom=="X"] <- 23
  tbl$chrom[tbl$chrom=="Y"] <- 24
  tbl$chrom <- as.numeric(tbl$chrom)
  tbl[order(tbl$chrom),]
}
getCentromeres <- function(genome="hg19"){
  mySession <- try(browserSession("UCSC"), silent=TRUE)
  # In case of failure, try another mirror
  if(inherits(mySession, "try-error"))
    mySession <- browserSession("UCSC",
                                url="http://genome-euro.ucsc.edu/cgi-bin/")
  genome(mySession) <- genome
  obj <- ucscTableQuery(mySession, table="gap")
  tbl <- getTable(obj)
  tbl <- tbl[tbl$type=="centromere", c("chrom", "chromStart", "chromEnd")]
  colnames(tbl)[2:3] <- c("centromerStart", "centromerEnd")
  .chrAsNum(tbl)
}
maskRegions <- function(GR, GR.mask){
  require(regioneR)
  GR <- GR[!GR %over% GR.mask]
  return(GR)
}
maskSubTelomers <- function(GR, genome){
  require(regioneR)
  GR <- GR[!GR %over% subTelomerEnd]
  return(GR)
}
