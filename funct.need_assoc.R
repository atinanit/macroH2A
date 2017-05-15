library(regioneR)
library(lattice)
library(colorRamps)
listRegionset<-function(cell.line,type,folder){
  a<-read.delim(paste(folder,"files.txt",sep=""))
  a<-as.character(a[,1])
  a<-a[grep(type,a)]
  a<-a[grep(cell.line,a)]
  a<-paste(folder,a,sep="")
  return(a)
}
listRegionsObject<-function(list){
  test<-list()
  for(i in 1:length(list)){
    download.file(list[i],"test.gz")
    test[i]<-toGRanges(read.delim(gzfile("test.gz")))
  }
  names(test)<-list
  return(test)
}
listOverlapPermTest<-function(Alist,Blist,ntimes,...){
  list.tabs<-list()
  list.pt<-list()
  for ( i in 1:length(Alist)){
    A<-Alist[[i]]
    new.names<-names(Blist)
    func.list <- createFunctionsList(FUN=numOverlaps, param.name="B", values=Blist)

    ptm <- proc.time()
    pt <- permTest(A=A, ntimes=ntimes, evaluate.function=func.list,
                   randomize.function=randomizeRegions,...)
    time<-proc.time() - ptm
    time<-time[3]/60
    print(paste0(" run in ",time,"  minute"))
    tab<-vector()
    for (j in 1:length(pt)){
      if(pt[[j]]$zscore>0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(Blist[[j]]))
        max.val<-length(Blist[[j]])
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(Blist[[j]]))
        zscore.std<-zscore.norm/max.z.norm
      }
      if(pt[[j]]$zscore<=0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(Blist[[j]]))
        max.val<-0
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(Blist[[j]]))
        zscore.std<-(-(zscore.norm/max.z.norm))
      }
      vec<-c(j,new.names[j],length(Blist[[j]]),pt[[j]]$zscore,zscore.norm,zscore.std,pt[[j]]$pval)
      tab<-rbind(vec,tab)
    }
    colnames(tab)<-c("order.id","name","nº region","z-score","norm.z-score","standard.z-score","p-value")
    print(tab)
    list.tabs[[i]]<-tab
    list.pt[[i]]<-pt
    names(list.tabs)[i]<-names(Alist)[i]
    names(list.pt)[i]<-names(Alist)[i]

  }
  return(list(list.tabs,list.pt))
}
# #randomizeRegions <- function(A, genome="hg19",  mask=NULL, non.overlapping=TRUE, per.chromosome=FALSE, ...) {
#
#   if(!hasArg(A)) stop("A is missing")
#   if(!is.logical(non.overlapping)) stop("non.overlapping must be logical")
#   if(!is.logical(per.chromosome)) stop("per.chromosome must be logical")
#
#
#   A <- toGRanges(A)
#
#   #The randomization of an empty region set, is an empty region set
#   if(length(A)<1) {
#     return(A)
#   }
#
#   gam <- getGenomeAndMask(genome=genome, mask=mask)
#   genome <- gam[["genome"]]
#   mask <- gam[["mask"]]
#
#   #use subtractRegions to get the masked genome
#   valid.regions <- toDataframe(subtractRegions(genome, mask))[,c(1,2,3)]
#   valid.regions[,1] <- as.character(valid.regions[,1])
#   names(valid.regions) <- c("chr", "start", "end")
#
#
#
#   if(per.chromosome) {
#     missing.chrs <- setdiff(levels(seqnames(A)), valid.regions[,1])
#     if(length(missing.chrs > 0)) {
#       stop(paste("The chromosomes **", paste(missing.chrs, sep=", "), "** of A are not present in the genome or are completely masked. It is not possible to use \"per.chromosome\" with this dataset.", sep=""))
#     }
#     #Split the regions in A and the valid regions per chromosomes and create the random regiosn separately per each chromosome
#     random.regions <- toGRanges(data.frame(chr=character(), start=numeric(), end=numeric(), stringsAsFactors=FALSE))
#     seqlevels(random.regions)<-seqlevels(A)
#     levels(seqnames(random.regions))<-seqlevels(A)
#
#
#     for(chr in seqlevels(A)) {
#       chr.A <- A[seqnames(A)==chr]
#       chr.valid <- valid.regions[valid.regions[,1]==chr,]
#       new.regions<-hybrid_randomizeRegions(chr.A, chr.valid, non.overlapping=non.overlapping)
#       seqlevels(new.regions)<-seqlevels(A)
#       levels(seqnames(new.regions))<-seqlevels(A)
#       random.regions <- c(random.regions, new.regions)
#     }
#   } else {
#     random.regions <- hybrid_randomizeRegions(A, valid.regions=valid.regions, non.overlapping=non.overlapping)
#   }
#
#   return(random.regions)
#
# }
# #hybrid_randomizeRegions <- function(A, valid.regions, non.overlapping=FALSE, max.retries=5) {
#
#   if(length(A)<1) {
#     return(A)
#   }
#
#   if(non.overlapping==FALSE) { #simply pass to the actual function
#     rr <- private_randomizeRegions(A=A, valid.regions=valid.regions, non.overlapping=non.overlapping, max.retries=max.retries)
#   } else {
#     #Create a set of regions allowing overlaps
#     rr <- private_randomizeRegions(A=A, valid.regions=valid.regions, non.overlapping=FALSE, max.retries=max.retries)
#
#
#     #Detect the overlaps between different regions
#     ov <- overlapRegions(A=rr, B=rr)
#     dups <- ov$type != "equal"   #TODO: take into account the situation where, by chance, two regions are exactly the same
#
#     if(length(which(dups))>0) {
#       pending <- toGRanges(ov[dups, c(1,2,3)]) #the regions taht overlapped something are still pending
#
#       toRemove <- which(paste0(seqnames(rr), "#",start(rr),"#", end(rr)) %in%  paste0(seqnames(pending), "#",start(pending),"#", end(pending)))
#
#       rr <- rr[-toRemove,] #the regions that overlapped nothing, are kept in the random regions set
#
#       #Add the already placed regions into the mask
#       valid.regions <- toDataframe(subtractRegions(valid.regions, rr))
#       #and place the remaining regions using the quadratic algorithm
#       #if the number of regions to place is high (>1000), call hybrid_randomizeRegions recursively, else, call the quadratic algorithm
#       #TODO: What if the regions do not fit into the genome? We should place a guard against taht situation limiting the number of recursive calls
#       if(length(pending)>500) {
#         rr2 <- hybrid_randomizeRegions(A=pending, valid.regions=valid.regions, non.overlapping=TRUE, max.retries=max.retries)
#       } else {
#         rr2 <- private_randomizeRegions(A=pending, valid.regions=valid.regions, non.overlapping=TRUE, max.retries=max.retries)
#       }
#       suppressWarnings(rr <- append(rr, rr2)) #supress warnings, since it will warn of potential different reference genome if we have just a few regions
#     }
#   }
#
#   return(rr)
#
# }
# #Curry <- function(FUN, ...) {
#   args <- match.call(expand.dots = FALSE)$...
#   args$... <- as.name("...")
#
#   env <- new.env(parent = parent.frame())
#
#   if (is.name(FUN)) {
#     fname <- FUN
#   } else if (is.character(FUN)) {
#     fname <- as.name(FUN)
#   } else if (is.function(FUN)){
#     fname <- as.name("FUN")
#     env$FUN <- FUN
#   } else {
#     stop("FUN not function or name of function")
#   }
#   curry_call <- as.call(c(list(fname), args))
#
#   f <- eval(call("function", as.pairlist(alist(... = )), curry_call))
#   environment(f) <- env
#   f
# }
# #createFunctionsList <- function(FUN, param.name, values, func.names=NULL) {
#   if(!hasArg(FUN)) stop("FUN is missing")
#   if(!hasArg(param.name)) stop("param.name is missing")
#   if(!is.character(param.name)) stop("param.name must be a character")
#   if(length(param.name)>1) stop("param.name must be a single character string")
#   if(!hasArg(values)) stop("values is missing")
#   if(!is.list(values)) stop("values must be a list")
#
#
#   if(is.null(func.names)) {
#     if(is.list(values) && !is.null(names(values))) {
#       func.names <- names(values)
#     } else {
#       func.names <- paste0("Function", c(1:length(values)))
#     }
#   }
#
#   curried.funcs <- list()
#   for(i in c(1:length(values))) {
#     curry.args <- list(FUN=FUN)
#     curry.args[[param.name]] <- values[[i]]
#     curried.funcs[[func.names[i]]] <- do.call("Curry", curry.args)
#   }
#   return(curried.funcs)
# }
# #private_randomizeRegions <- function(A, valid.regions, non.overlapping=FALSE, max.retries=5) {
#
#   if(length(A)<1) {
#     return(A)
#   }
#
#   #Convenience function to map back the position in the vector of valid positions to the original regions
#   map.back <- function(pos) {
#     num.region <- which(len.valid >= pos & len.valid != 0)[1]  #We can be sure al least one len.valid will be >= than pos
#     if(num.region>1) pos <- pos - len.valid[num.region-1]
#     #start <- data.frame(chr=valid.regions[num.region, 1], start=valid.regions[num.region, 2]+pos)
#     #new.pos <- data.frame(chr=as.character(seqnames(valid.regions)[num.region]), start=start(valid.regions)[num.region]+pos)
#     new.pos <- list(num.region=num.region, chr=as.character(valid.regions[num.region, 1]), start=valid.regions[num.region, 2]+pos)
#     return(new.pos)
#   }
#
#
#   #TODO: Should we really retry? or just giveup immediately?
#   #Set up a retry system in case the region set does not fit into the genome with a given configuration
#   original.valid.regions <- valid.regions
#   for(ntry in 1:max.retries) {
#
#     #reset the valid.regions
#     valid.regions <- original.valid.regions
#
#     #TODO: Sort from largest to smallest to improve the probability of finding a place for all regions in highly fragmented genomes
#
#     failed <- FALSE
#     first <- TRUE
#     #Bernat
#     new.chr <- new.start <- new.end <- c()
#     for(i in 1:length(A)) {
#
#       len <- width(A)[i]
#       #CHECK: a -1 is not necessary when using the dataframe fr valid.regions
#       len.valid <- valid.regions[,3]-valid.regions[,2] - len -1 #The region cannot start in the last len positions of the valid region. It would not fit
#       #len.valid <- width(valid.regions) - len - 1 #The region cannot start in the last len positions of the valid region. It would not fit
#       len.valid[len.valid<0] <- 0
#       if(length(len.valid)>1) {
#         for(j in 2:length(len.valid)) {
#           len.valid[j] <- len.valid[j] + len.valid[j-1]
#           #print(paste(j,len.valid[j]))
#         }
#       }
#
#       if(max(len.valid)==0) { #If theres no space left to put this region, stop this try
#         failed = TRUE
#         break
#       }
#
#       rand.num <- round(runif(1)*len.valid[length(len.valid)])
#
#       new.pos <- map.back(rand.num)
#
#
#       #New without dataframe
#       new.chr <- c(new.chr, new.pos[["chr"]])
#       new.start <- c(new.start, new.pos[["start"]])
#       new.end <- c(new.end, new.pos[["start"]] + len)
#
#       #Finally, if non overlapping regions are required, remove the last region from the valid regions set
#       #Note: the non overlapping processing changes the order of the regions, but that does not have any impact on the randomness
#       if(non.overlapping) {
#         #Substitute the current region (old, num.region) by two new ones, substracting the new random region just created
#
#         old.end <- valid.regions[new.pos[["num.region"]],]
#         old.end[1,"start"] <- new.pos[["start"]]+len+1
#
#         valid.regions[new.pos[["num.region"]], "end"] <- new.pos[["start"]]-1
#
#         valid.regions <- rbind(valid.regions, old.end)
#
#         #If the newly created region was at the exact start position or end position of a valid region, this will create a valid region with size -1. Remove them
#         #TODO: we could check this condition before instead of creating and removing.
#         negative.lengths <- (valid.regions[,"end"]-valid.regions[,"start"])<0
#         if(any(negative.lengths)) {
#           valid.regions <- valid.regions[-negative.lengths,]
#         }
#
#       }
#
#     }
#
#     if(!failed) {
#       random.regions <- data.frame(chr=new.chr, start=new.start, end=new.end)
#       return(toGRanges(random.regions))
#     } else {
#       warning("It was not possible to create the random region set because there was no space available. Retrying.")
#     }
#   }
#   #If are here, all the retries have failed to produce a valid region set. Raise an error and stop.
#   stop(paste("It was not possible to create a valid random region set after ", max.retries, " retries. This might be due to the regions in A covering most of the available genome. Allowing overlapping random regions could solve this problem.", sep=""))
# }
clustMatrix<-function(mat,...){
  d<-dist(mat)
  hc<-hclust(d,...)
  ind<-hc$order

  d1<-dist(t(mat))
  hc1<-hclust(d1,...)
  ind1<-hc1$order

  mat1<-mat[ind,ind1]
  return(mat1)
}
clustMatrix.double<-function(mat,...){

  ordx<-order(colnames(mat))
  ordy<-order(rownames(mat))
  mat<-mat[ordy,ordx]
  d<-dist(mat)
  hc<-hclust(d,...)
  ind<-hc$order

  mat1<-mat[ind,ind]
  return(mat1)
}
listOverlapPermTest.true<-function(Alist,Blist,ntimes,...){
  list.tabs<-list()
  list.pt<-list()
  for ( i in 1:length(Alist)){
    A<-Alist[[i]]
    new.names<-names(Blist)
    func.list <- createFunctionsList(FUN=numOverlaps, param.name="B", values=Blist)

    ptm <- proc.time()
    pt <- permTest(A=A, ntimes=ntimes, evaluate.function=func.list,
                   randomize.function=randomizeRegions,count.once=T,per.chromosome=T)
    time<-proc.time() - ptm
    time<-time[3]/60
    print(paste0(" run in ",time,"  minute"))
    tab<-vector()
    for (j in 1:length(pt)){
      if(pt[[j]]$zscore>0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
        max.val<-min(length(A),length(Blist[[j]]))
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(A))
        zscore.std<-zscore.norm/max.z.norm
      }
      if(pt[[j]]$zscore<=0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
        max.val<-0
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(A))
        zscore.std<-(-(zscore.norm/max.z.norm))
      }
      vec<-c(j,new.names[j],length(Blist[[j]]),pt[[j]]$zscore,zscore.norm,zscore.std,pt[[j]]$pval,pt[[j]]$observed,
             mean(pt[[j]]$permuted), sd(pt[[j]]$permuted))
      tab<-rbind(vec,tab)
    }
    colnames(tab)<-c("order.id","name","nº region","z-score","norm.z-score","standard.z-score","p-value","n.overlaps","mean.perm","sd.perm")
    print(tab)
    list.tabs[[i]]<-tab
    #list.pt[[i]]<-pt
    names(list.tabs)[i]<-names(Alist)[i]
    #names(list.pt)[i]<-names(Alist)[i]

  }
  return(list.tabs)
}
matListOver<-function(listOverlapPermTest.obj,zs.type="std"){
  A.obj<-listOverlapPermTest.obj
  mat<-vector()
  for (i in 1:length(A.obj)){
    if (zs.type=="std"){
      mat<-cbind(mat,as.numeric(A.obj[[i]][,6]))
    }
    if (zs.type=="norm"){
      mat<-cbind(mat,as.numeric(A.obj[[i]][,5]))
    }
  }
  colnames(mat)<-names(A.obj)
  rownames(mat)<-A.obj[[1]][,2]
  mat<-as.matrix(mat)
  return(mat)
}
permTest <- function(A, ntimes=100, randomize.function, evaluate.function, alternative="auto", min.parallel=1000, force.parallel=NULL, randomize.function.name=NULL, evaluate.function.name=NULL, verbose=FALSE, ...) {

  #check arguments
  alternative<-match.arg(alternative,c("less","greater", "auto"))
  if(!hasArg(A)) stop("A is missing")
  if(!is.numeric(ntimes)) stop("ntime must be numeric")
  if(!hasArg(randomize.function)) stop("randomize.function is missing")
  if(!is.function(randomize.function)) stop("randomize.function must be a function")
  if(!hasArg(evaluate.function)) stop("evaluate.function is missing")
  if(!(is.function(evaluate.function) | is.list(evaluate.function))) stop("evaluate.function must be a function")
  if(!is.numeric(min.parallel)) stop("min.parallel must be numeric")
  if(ntimes<100) print(paste0("Note: The minimum p-value with only ",ntimes," permutations is ",1/(ntimes+1),". You should consider increasing the number of permutations."))

  A <- toGRanges(A) #does nothing if already a GRanges object


  if(!is.null(force.parallel)) {
    doParallel <- force.parallel
  } else {
    doParallel <- (length(A)*ntimes > min.parallel)
  }


  #Evaluation Function: get the function name and convert to list if its not yet
  if(!is.list(evaluate.function)) { #if it's a single function
    if(is.null(evaluate.function.name)) {
      evaluate.function.name <- as.character(match.call()["evaluate.function"])
    }
    ef <- list()
    ef[[evaluate.function.name]] <-  evaluate.function
    evaluate.function <- ef
  } else { #if it's a list of functions
    if(!is.null(evaluate.function.name)) { #if names were explicitely provided
      names(evaluate.function) <- evaluate.function.name
    } else { #try to leave the current names or create new ones if no names present
      if(is.null(names(evaluate.function))) {
        names(evaluate.function) <- paste0("Function", c(1:length(evaluate.function)))
      }
    }
  }

  #Randomization Function: Get a name
  if(is.null(randomize.function.name)) {
    randomize.function.name <- match.call()["randomize.function"]
  }

  #Start the permutation test
  #compute the evaluation function(s) using the original region set A
  original.evaluate <- sapply(c(1:length(evaluate.function)), function(i,...) {return(evaluate.function[[i]](A,...))}, ...)

  if(!is.numeric(original.evaluate)) {
    stop(paste0("The evaluation function must return a numeric value but it returned an object of class ", class(original.evaluate)))
  }

  if(verbose) {
    #WARNING: to give some visual information about the computation done, we use a Progress Bar. However, since the GUI will not be updated
    #whilst in multi core processing (mclapply), we partition the work in chunks and repeatedly call the parallel computation, once per chunk.
    #Between chunks, the progess bar will be updated. The problem is that the total computation time might increase due to waiting and joining
    #multiple times.
    #Create the progress bar

    pb <- txtProgressBar(min = 0, max = ntimes, style = 3)
    setTxtProgressBar(pb, 0)
  }

  #define the function to create and evaluate the random sets
  randomize_and_evaluate <- function(foo, ...) {
    #randomize
    randomA <- randomize.function(A,...)
    #evaluate the random region set
    if(verbose) {
      setTxtProgressBar(pb, foo)
    }

    #compute the evaluation function(s) using the RANSOMIZED region set randomA
    rand.evaluate <- sapply(c(1:length(evaluate.function)), function(i, ...) {return(evaluate.function[[i]](randomA,...))}, ...)

    return(rand.evaluate)
  }

  #create the random sets and evaluate them
  if(doParallel) {
    if(verbose) { #if verbose, we will do the computations in chunks and update the progress bar in between
      random.evaluate <- numeric()
      chunk.size <- max(round(ntimes/100+1), 10)
      e <- 0
      done <- FALSE
      while(!done) {
        s <- e + 1
        e <- s + chunk.size
        if(e >= ntimes) {
          e <- ntimes
          done <- TRUE
        }
        random.evaluate <- c(random.evaluate, do.call(rbind, mclapply(c(s:e), randomize_and_evaluate, ...)))
        setTxtProgressBar(pb, e)
      }
    } else { #if not verbose, just do it
      random.evaluate <- do.call(rbind, mclapply(c(1:ntimes), randomize_and_evaluate, ...))
    }
  } else {
    random.evaluate <- do.call(rbind, lapply(c(1:ntimes), randomize_and_evaluate, ...))
  }


  #The simulation process has finished. Now build a permTestResults object for each evaluate.function
  results <- list()
  for(i in c(1:length(evaluate.function))) {
    #Get the data for the i-th function
    func.name <- names(evaluate.function)[i]
    orig.ev <- original.evaluate[i]
    rand.ev <- random.evaluate[,i]


    #warn if any NA
    num.nas <- length(which(is.na(rand.ev)))
    if(num.nas > 0) warning(paste0(num.nas, " iterations returned NA's. Only ", ntimes-num.nas ," iterations have been used to compute the p-value."))

    #decide the alternative if alternative == "auto"
    if(alternative == "auto") {
      alt <- ifelse(orig.ev < mean(rand.ev, na.rm=TRUE), "less", "greater")
    } else {
      alt <- alternative
    }

    #Compute the p-value
    if (alt == "less") {
      pval <- (sum(orig.ev > rand.ev, na.rm=TRUE) + 1) / (ntimes - num.nas + 1)
    } else { #alt == "greater"
      pval <- (sum(orig.ev < rand.ev, na.rm=TRUE) + 1) / (ntimes - num.nas + 1)
    }
    #if the original alternative was not the best one, suggest the user to change it
    if(alternative=="greater" & orig.ev<mean(rand.ev,na.rm=TRUE)) message("Alternative is greater and the observed statistic is less than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
    if(alternative=="less" & orig.ev>mean(rand.ev,na.rm=TRUE)) message("Alternative is less and the observed statistic is greater than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")

    #Compute the z-score
    if(orig.ev == 0 & all(rand.ev == 0)){ #If everything is 0, warning and "empty" results
      warning(paste0("All permuted values and the original evaluation value are equal to 0. Z-score cannot be computed."))
      pval <- 1
      zscore <- NA
    } else{
      zscore <- round((orig.ev - mean(rand.ev, na.rm=TRUE)) / sd(rand.ev, na.rm=TRUE), 4)
    }
    if(!is.finite(zscore)){ #if all evaluations are equal, the sd is 0 and the z-score is infinite
      warning(paste0("All permuted values are equal to ", rand.ev[1], ". Z-score is infinite."))
    }

    #Create the permTestResults object
    res<-list(pval=pval, ntimes=ntimes, alternative=alt, observed=orig.ev, permuted=rand.ev, zscore=zscore,
              evaluate.function=evaluate.function[[i]], evaluate.function.name=func.name,
              randomize.function=randomize.function, randomize.function.name=randomize.function.name)
    class(res) <- "permTestResults"
    results[[func.name]] <- res
  }


  return(results)

}
listOverlapPermTest.true.circ<-function(Alist,Blist,ntimes,...){
  list.tabs<-list()
  list.pt<-list()
  for ( i in 1:length(Alist)){
    A<-Alist[[i]]
    A<-filterChromosomes(A)
    new.names<-names(Blist)
    func.list <- createFunctionsList(FUN=numOverlaps, param.name="B", values=Blist)

    ptm <- proc.time()
    pt <- permTest(A=A, ntimes=ntimes, evaluate.function=func.list,
                   randomize.function=circularRandomizeRegions,...)
    time<-proc.time() - ptm
    time<-time[3]/60
    print(paste0(" run in ",time,"  minute"))
    tab<-vector()
    for (j in 1:length(pt)){
      if(pt[[j]]$zscore>0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
        max.val<-min(length(A),length(Blist[[j]]))
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(A))
        zscore.std<-zscore.norm/max.z.norm
      }
      if(pt[[j]]$zscore<=0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
        max.val<-0
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(A))
        zscore.std<-(-(zscore.norm/max.z.norm))
      }
      vec<-c(j,new.names[j],length(Blist[[j]]),pt[[j]]$zscore,zscore.norm,zscore.std,pt[[j]]$pval,pt[[j]]$observed,
             mean(pt[[j]]$permuted), sd(pt[[j]]$permuted))
      tab<-rbind(vec,tab)
    }
    colnames(tab)<-c("order.id","name","nº region","z-score","norm.z-score","standard.z-score","p-value","n.overlaps","mean.perm","sd.perm")
    print(tab)
    list.tabs[[i]]<-tab
    #list.pt[[i]]<-pt
    names(list.tabs)[i]<-names(Alist)[i]
    #names(list.pt)[i]<-names(Alist)[i]

  }
  return(list.tabs)
}
ReduxAndCompareMatrix<-function(mat.A,mat.B,nameA="",nameB="",...){

  colnames(mat.A)<-names.A<-sub(".*-", "", colnames(mat.A))
  colnames(mat.B)<-names.B<-sub(".*-", "", colnames(mat.B))
  rownames(mat.A)<-sub(".*-", "",rownames(mat.A))
  rownames(mat.B)<-sub(".*-", "",rownames(mat.B))

  mat1<-mat.A[which(rownames(mat.A)%in%rownames(mat.B)),which(colnames(mat.A)%in%colnames(mat.B))]
  mat1<-mat1[!duplicated(rownames(mat1)),!duplicated(colnames(mat1))]
  mat2<-mat.B[which(rownames(mat.B)%in%rownames(mat.A)),which(colnames(mat.B)%in%colnames(mat.A))]
  mat2<-mat2[!duplicated(rownames(mat2)),!duplicated(colnames(mat2))]
  mat1<-mat1[order(rownames(mat1)),order(colnames(mat1))]
  mat2<-mat2[order(rownames(mat2)),order(colnames(mat2))]

  mat1<-clustMatrix(mat1)
  mat2<-mat2[match(rownames(mat1),rownames(mat2)),match(colnames(mat1),colnames(mat2))]
  mat3<-mat1-mat2

  mat.tot<-list(mat1,mat2,mat3)
  names(mat.tot)<-c(nameA,nameB,paste0(nameA,"-",nameB))
  return(mat.tot)

}
TabToListHMM<-function(tab,vect.test){
  lis<-list()
  names<-as.vector(unique(vect.test))
  for(i in 1:length(names)){
    lis[[i]]<-toGRanges(tab[which(vect.test==names[i]),])
  }
  names(lis)<-names

  for(i in 1:length(lis)){
    if (length(lis[[i]])>50000){
      lis[[i]]<-lis[[i]][sample(length(lis[[i]]),50000)]
    }
  }

  for(i in 1:length(lis)){
    names(lis)[i]<-paste0(i,"-",names(lis)[i])
  }
  return(lis)
}
listOverlapPermTest.resampling<-function(Alist,Blist,ntimes,universe="auto",mc.cores=4,...){
  list.tabs<-list()
  list.pt<-list()
  if (universe=="auto"){
    DFF<-data.frame()
    for(i in 1:length(Alist)){
      DF<-toDataframe(Alist[[i]])[,1:3]
      DFF<-rbind(DFF,DF)
    }
    universe<-DFF
  }
  for ( i in 1:length(Alist)){
    A<-Alist[[i]]
    new.names<-names(Blist)
    func.list <- createFunctionsList(FUN=numOverlaps, param.name="B", values=Blist)

    ptm <- proc.time()
    pt <- permTest(A=A, ntimes=ntimes, evaluate.function=func.list, count.once=T,
                   randomize.function=resampleRegions,universe=universe,per.chromosome=F,mc.cores=mc.cores)
    time<-proc.time() - ptm
    time<-time[3]/60
    print(paste0(" run in ",time,"  minute"))
    tab<-vector()
    for (j in 1:length(pt)){
      if(pt[[j]]$zscore>0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
        max.val<-min(length(A),length(Blist[[j]]))
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(A))
        zscore.std<-zscore.norm/max.z.norm
      }
      if(pt[[j]]$zscore<=0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
        max.val<-0
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(A))
        zscore.std<-(-(zscore.norm/max.z.norm))
      }
      vec<-c(j,new.names[j],length(Blist[[j]]),pt[[j]]$zscore,zscore.norm,zscore.std,pt[[j]]$pval,pt[[j]]$observed,
             mean(pt[[j]]$permuted), sd(pt[[j]]$permuted))
      tab<-rbind(vec,tab)
    }


    tab<-cbind(tab,p.adjust(as.numeric(tab[,7]),method="BH"))
    colnames(tab)<-c("order.id","name","nº region","z-score","norm.z-score","standard.z-score","p-value","n.overlaps","mean.perm","sd.perm","adj.p-value")
    print(tab)
    list.tabs[[i]]<-tab
    #list.pt[[i]]<-pt
    names(list.tabs)[i]<-names(Alist)[i]
    #names(list.pt)[i]<-names(Alist)[i]

  }
  return(list.tabs)
}
plotMatrix<-function(mat,palette="bwr",type="normal",name="",limit=TRUE,gr.range=c(-1.2,1.2)){

  atx<-seq(min(mat),max(mat),length.out=100)

  if(limit==TRUE){
    atx<-seq(gr.range[1],gr.range[2],length.out=100)
  } else {

  }
  if (palette=="bwr"){rgb.palette <- colorRampPalette(c("darkblue","blue","white","red", "darkred"),space = "rgb")}
  if (palette=="gbr"){rgb.palette <- colorRampPalette(c("green","black","red"),space = "rgb")}
  if (palette=="bwy"){rgb.palette <- colorRampPalette(c("black","white","yellow"),space = "rgb")
  }
  trellis.par.set(regions=list(col=rgb.palette(100)))
  if (type=="normal"){
    pl<-levelplot(mat,at=atx, scale=list(x=list(rot=45)),main=name, aspect="fill")
  }
  if (type=="compare"){
    pl<-levelplot(mat, scale=list(x=list(rot=45)),main=name, aspect="fill")
  }
  return(pl)
}
clustMatrix.double.2<-function(mat1,method="complete",...){
  mat<-mat1
  mat<-mat^1
  ordx<-order(colnames(mat))
  ordy<-order(rownames(mat))
  mat<-mat[ordy,ordx]
  mat1<-mat1[ordy,ordx]
  d<-dist(mat)
  hc<-hclust(d,method=method)
  ind<-hc$order

  mat1<-mat1[ind,ind]
  return(mat1)
}
selectCalls <- function(calls){
  for (i in 1:length(calls)){
    if (length(calls[[i]]) > 5000){
      ind <- sample(length(calls[[i]]), 5000)
      calls[[i]] <- calls[[i]][ind]
    }
  }
  return (calls)
}
globalAssociations <- function(call, GRanges, count.once = TRUE, lz = TRUE, ntimes = 100, mc.cores=8, window = 10000, step = 100, genome = "hg19", evaluate.func = "numOverlaps"){
  if (evaluate.func == "numOverlaps"){
    function.list <- createFunctionsList(FUN = numOverlaps, param.name = "B", values = GRanges)
    m_rmsk.pt <- permTest(A = call,
                          ntimes = ntimes,
                          evaluate.function = function.list,
                          randomize.function = randomizeRegions,
                          genome = genome,
                          per.chromosome = TRUE,
                          count.once = count.once,
                          mc.cores = mc.cores)
  }
  if (evaluate.func == "ATGC.eval"){
    m_rmsk.pt <- permTest(A = call,
                          ntimes = ntimes,
                          evaluate.function = ATGC.eval,
                          randomize.function = randomizeRegions,
                          genome = genome,
                          per.chromosome = TRUE,
                          count.once = count.once,
                          mc.cores = mc.cores)
  }
  if (lz == TRUE){
    lz_m_rmsk <- lzglobalAs(call, m_rmsk.pt, window = window, step = step, count.once=count.once)
    return (list(m_rmsk.pt,lz_m_rmsk))
  } else {
    return (m_rmsk.pt)
  }
}
countBase<-function(base,string){
  p <- base
  s2 <- gsub(p,"",string)
  numOcc <- nchar(string) - nchar(s2)
  return(numOcc)
}
perBase<-function(string){
  tot<-nchar(string)
  nA<-countBase("A",string)
  nT<-countBase("T",string)
  nC<-countBase("C",string)
  nG<-countBase("G",string)
  pA<-nA/tot
  pT<-nT/tot
  pC<-nC/tot
  pG<-nG/tot
  res<-data.frame(pA+pT,pC+pG)
  return(res)
}
parse.bed<-function(bed,chr,st=0,en=1e20){
  parsed<-bed[bed[,1]==chr,]
  parsed<-parsed[parsed[,2]>=st,]
  parsed<-parsed[parsed[,3]<=en,]
  return(parsed)
}
ATGC.eval<-function(A,type="AT",genome="hg19",...){
  A<-toDataframe(A)
  ATGCper<-data.frame()
  for(i in nrow(A)){
    if (genome=="hg19"){gs<-toString(getSeq(Hsapiens, as.character(A[i,1]), as.integer(A[i,2]), as.integer(A[i,3])))}
    if (genome=="mm9"){gs<-toString(getSeq(Mmusculus, as.character(A[i,1]), as.integer(A[i,2]), as.integer(A[i,3])))}
    ATGCper<-rbind(ATGCper,perBase(gs))
  }
  if (type=="GC"){return(as.numeric(colMeans(ATGCper)[1]))}
  if (type=="AT"){return(as.numeric(colMeans(ATGCper)[2]))}
}
lzglobalAs <- function(call, m_rmsk.pt, window = 10000, step = 100, count.once = TRUE){
  lz_m_rmsk <- list()
  for(i in 1:length(m_rmsk.pt)){
    lz_m_rmsk[[i]] <- localZScore(A = call,  pt = m_rmsk.pt[[i]], window = window, step = step, count.once = count.once)
    names(lz_m_rmsk)[i]<-names(m_rmsk.pt[i])
  }
  return(lz_m_rmsk)
}
matCreate <- function(call, lz_m_rmsk){
  lz_m_rmsk_shifted <- vector()
  lz.tabs <- lz_m_rmsk[2][[1]]
  for(i in 1:length(lz.tabs)){
    lz_m_rmsk_shifted <- rbind(lz_m_rmsk_shifted, lz.tabs[[i]]$shifted.z.scores/sqrt(length(call)))
  }
  #lz_m_rmsk_shifted[(lz_m_rmsk_shifted < cutoff) & (lz_m_rmsk_shifted > -cutoff)] <- 0
  rownames(lz_m_rmsk_shifted) <- names(lz.tabs)
  lz_m_rmsk_shifted <- as.matrix(lz_m_rmsk_shifted)
  return (lz_m_rmsk_shifted)
}
matCreate_matrix <- function(call, lz_m_rmsk){
  # argument is a matrix
  lz_m_rmsk_shifted <- vector()
  #lz.tabs <- lz_m_rmsk[2][[1]]
  for(i in 1:length(lz_m_rmsk)){
    lz_m_rmsk_shifted <- rbind(lz_m_rmsk_shifted, lz_m_rmsk[[i]]$shifted.z.scores/sqrt(length(call)))
  }
  rownames(lz_m_rmsk_shifted) <- names(lz_m_rmsk)
  lz_m_rmsk_shifted <- as.matrix(lz_m_rmsk_shifted)
  return (lz_m_rmsk_shifted)
}
clustMatrix.column<-function(mat,...){
  d<-dist(mat)
  hc<-hclust(d)
  ind<-hc$order

  mat1<-mat[ind,]
  return(list(mat1, ind))
}
getMatrixAssociation <- function(call, lz_m_rmsk, name = "", mat.plot = FALSE, ind = NULL, limit=TRUE, gr.range=c(-0.8,0.8)){
  lz_m_rmsk_shifted <- matCreate(call, lz_m_rmsk)
  if(is.null(ind)){
    list_lz_m_rmsk_shifted.HC_ind <- clustMatrix.column(lz_m_rmsk_shifted)
    lz_m_rmsk_shifted.HC <- list_lz_m_rmsk_shifted.HC_ind[[1]]
  } else {
    lz_m_rmsk_shifted.HC <- lz_m_rmsk_shifted[ind,]
  }
  names <- rownames(lz_m_rmsk_shifted.HC)
  if(mat.plot == TRUE){
    matplot(t(lz_m_rmsk_shifted.HC), type="l", col = primary.colors(16), main = name)
    #legend("topleft",legend=c(names), col = primary.colors(16), lty=1:16, cex=0.8)
  }
  plot.matrix <- plotMatrix(t(lz_m_rmsk_shifted.HC),palette="gbr",type="normal",name = name,limit=TRUE,gr.range=gr.range)
  if(is.null(ind)){
    return (list(list_lz_m_rmsk_shifted.HC_ind, lz_m_rmsk_shifted.HC, plot.matrix))
  } else {
    return (list(lz_m_rmsk_shifted.HC, plot.matrix))
  }
}
getMatrixAssociation_matrix <- function(call, lz_m_rmsk, name = "", mat.plot = FALSE, ind = NULL, limit=TRUE, gr.range=c(-0.8,0.8)){
  lz_m_rmsk_shifted <- matCreate_matrix(call, lz_m_rmsk)
  if(is.null(ind)){
    list_lz_m_rmsk_shifted.HC_ind <- clustMatrix.column(lz_m_rmsk_shifted)
    lz_m_rmsk_shifted.HC <- list_lz_m_rmsk_shifted.HC_ind[[1]]
  } else {
    lz_m_rmsk_shifted.HC <- lz_m_rmsk_shifted[ind,]
  }
  names <- rownames(lz_m_rmsk_shifted.HC)
  if(mat.plot == TRUE){
    matplot(t(lz_m_rmsk_shifted.HC), type="l", col = primary.colors(16), main = name)
    #legend("topleft",legend=c(names), col = primary.colors(16), lty=1:16, cex=0.8)
  }
  plot.matrix <- plotMatrix(t(lz_m_rmsk_shifted.HC),palette="gbr",type="normal",name = name,limit=TRUE,gr.range=gr.range)
  if(is.null(ind)){
    return (list(list_lz_m_rmsk_shifted.HC_ind, lz_m_rmsk_shifted.HC, plot.matrix))
  } else {
    return (list(lz_m_rmsk_shifted.HC, plot.matrix))
  }
}
changenames <- function(file){
  for(i in 1:length(file)){
    names(file)[i]<-strsplit(names(file)[i],"UniPk")[[1]][1]
    names(file)[i]<-paste0(i,"-",strsplit(names(file)[i],"Hepg2")[[1]][2])
  }
  return (file)
}
chromAutosomal <- function(list, organism = "hg"){
  for (i in 1:length(list)){
    list[[i]] <- filterChromosomes(list[[i]], organism = organism, chr.type = "autosomal")
  }
  return (list)
}
filterFunction <- function(listGRanges, vecmin = 500){
  vec <- vector()
  for(i in 1:length(listGRanges)){
    vec[i] <- length(listGRanges[[i]])
  }
  listGRanges <- listGRanges[vec > vecmin]
  return(listGRanges)
}
dirtInputCleaner<-function(test,test.q="0.995",input.call){
  call.test<-gscreen(paste0(test," > ",test.q),iterator=75)
  call.test.input<-extendRegions(input.call,extend.start=75,extend.end=75)
  ind<-(which(overlapRegions(call.test,call.test.input,only.boolean=T)))
  call.test<-toGRanges(call.test[-ind,])
  return(call.test)
}
multiple_merges <- function(Region){
  merges_Region <- Region[[1]]
  for(i in 2:length(Region) - 1){
    test <- mergeRegions(Region[[i]], Region[[i+1]])
    merges_Region <- mergeRegions(merges_Region,test)
  }
  return(merges_Region)
}
