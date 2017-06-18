library(regioneR)
library(dplyr)
library(rtracklayer)
options(scipen = 999)

'%!in%' <- function(x,y)!('%in%'(x,y))

getMaxScoreRegion <- function(get.score.regions.file, GR.row){
  #Need max scoring of a region
  require(rtracklayer)
  bedGraph <- import(get.score.regions.file, which = GRanges(GR.row[1],
                                                             IRanges(
                                                               as.numeric(
                                                                 GR.row[2]),
                                                               as.numeric(
                                                                 GR.row[3]))))
  return(max(bedGraph$score))
}
getMaxScoreRegionBed <- function(bed, get.score.regions.file){
  #Need scoring of a BED file
  require(dplyr)
  bedGraph <- bed %>% rowwise() %>% mutate(score = getMaxScoreRegion(
    get.score.regions.file, chr, start, end))
  return(GRanges(data.frame(bedGraph)))
}
reBin<-function(GR, bin=0){
  fst<-start(GR)[1]
  lst<-end(GR)[length(GR)]
  rebin<-seq(fst,lst,bin)

  bed.rebin<-data.frame(as.character(seqnames(GR)[1]),rebin,(rebin+bin)-1)
  score<-vector()
  for( i in 1:nrow(bed.rebin)){
    ov<-overlapRegions(bed.rebin[i,],GR,colB=1)
    score[i]<-mean(ov$score)
    if(nrow(ov)==0){score[i]<-0}
  }

  bed.rebin$score<-score
  return(toGRanges(bed.rebin))
}
modelingScore <- function(bedGraph, bin = 0, spar = 0.4){
  if(bin != 0){
    bedGraph <- reBin(bedGraph,bin=bin)
  }
  if (length(bedGraph) >= 4){
    a <- smooth.spline(bedGraph$score,spar=spar)
    mod <- a$fit$coef/max(a$fit$coef)
  } else {
    mod <- 0
  }
  return(mod)
}

getScoringRegionIR <- function(get.score.regions.file, GR.row, modA.reference,
                               bin = 0, spar = 0.4, reverse = FALSE){
  require(rtracklayer)
  bedGraph <- import(get.score.regions.file, which = GRanges(GR.row[1],
                                                             IRanges(
                                                               as.numeric(
                                                                 GR.row[2]),
                                                               as.numeric(
                                                                 GR.row[3]))))

  modB <- modelingScore(bedGraph, bin = bin, spar = spar)
  ext <- min(length(modA.reference), length(modB))
  modA.reference <- modA.reference[1:ext]
  modB <- modB[1:ext]
  if (reverse == TRUE){
    cor.forward <- cor(modA.reference, modB, method="pearson")
    cor.rev <- cor(rev(modA.reference), modB, method="pearson")
    correlation <- list(cor.forward, cor.rev)
    return(correlation)
  } else {
    cor.forward <- cor(modA.reference, modB, method="pearson")
    return(cor.forward)
  }
}
getScoringRegion <- function(bedGraph, modA.reference,
                             bin = 0, spar = 0.4, reverse = FALSE){

  if (reverse == TRUE){
    correlation <- correlationMod(modA.reference, modB, method.cor = "pearson",
                                  reverse = TRUE)
    return(correlation)
  } else {
    cor.forward <- correlationMod(modA.reference, modB, method.cor = "pearson",
                                  reverse = FALSE)
    return(cor.forward)
  }
}
correlationMod <- function(modA.reference, modB, method.cor = "pearson",
                           reverse = FALSE){
  ext <- min(length(modA.reference), length(modB))
  modA.reference <- modA.reference[1:ext]
  modB <- modB[1:ext]
  if (reverse == TRUE){
    cor.forward <- cor(modA.reference, modB, method=method.cor)
    cor.rev <- cor(rev(modA.reference), modB, method=method.cor)
    return(list(cor.forward, cor.rev))
  } else {

  }
}
correlationModReverse <- function(modA.reference, modB,
                                  method.cor = "pearson"){
  ext <- min(length(modA.reference), length(modB))
  modA.reference <- rev(modA.reference[1:ext])
  modB <- modB[1:ext]

  cor.reverse <- cor(modA.reference, modB, method=method.cor)
  return(cor.reverse)
}
correlationIRanges <- function(get.score.regions.file, GR.row, modA.reference,
                               bin = 0, spar = 0.4, reverse = FALSE){
  require(regioneR)
  modB <- getScoringRegion(get.score.regions.file = get.score.regions.file,
                           GR.row = GR.row, bin = bin, spar = spar)
  if (reverse == TRUE){
    correlation <- correlationMod(modA.reference, modB, method.cor = "pearson",
                                  reverse = TRUE)
    return(correlation)
  } else {
    cor.forward <- correlationMod(modA.reference, modB, method.cor = "pearson",
                                  reverse = FALSE)
    return(cor.forward)
  }
}
associatedCorrelation <- function(get.score.regions.file, GR, modA.reference,
                                  bin = 0, spar = 0.4, get.max.score = FALSE,
                                  reverse = FALSE){
  require(dplyr)
  #require(parallel)
  #cl <- detectCores() - 2
  #cluster <- makeCluster(cl)
  if (reverse == FALSE){
    if (get.max.score == TRUE){
      correlation <- apply(GR, 1, getScoringRegionIR,
                           get.score.regions.file = get.score.regions.file,
                           modA.reference = modA.reference,
                           bin = bin,
                           spar = spar,
                           reverse = FALSE)
      GR$cor.forward <- correlation
      GR$max.score <- apply(GR, 1, getMaxScoreRegion,
                            get.score.regions.file = get.score.regions.file)
    } else {
      correlation <- apply(GR, 1, getScoringRegionIR,
                           get.score.regions.file = get.score.regions.file,
                           modA.reference = modA.reference,
                           bin = bin,
                           spar = spar,
                           reverse = FALSE)
      GR$cor.forward <- correlation
    }
  } else { #reverse TRUE
    if (get.max.score == TRUE){
      correlation <- apply(GR, 1, getScoringRegionIR,
                           get.score.regions.file = get.score.regions.file,
                           modA.reference = modA.reference,
                           bin = bin,
                           spar = spar,
                           reverse = TRUE)
      GR$cor.forward <- unlist(correlation)[c(TRUE, FALSE)]
      GR$cor.reverse <- unlist(correlation)[c(FALSE, TRUE)]
      GR$max.score <- apply(GR, 1, getMaxScoreRegion,
                            get.score.regions.file = get.score.regions.file)

    } else {
      correlation <- apply(GR, 1, getScoringRegionIR,
                           get.score.regions.file = get.score.regions.file,
                           modA.reference = modA.reference,
                           bin = bin,
                           spar = spar,
                           reverse = TRUE)
      GR$cor.forward <- unlist(correlation)[c(TRUE, FALSE)]
      GR$cor.reverse <- unlist(correlation)[c(FALSE, TRUE)]
    }
  }
  #stopCluster(cluster)
  return (GR)
}
