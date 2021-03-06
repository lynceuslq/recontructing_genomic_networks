#!/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Arguments must be supplied", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  
  print(args)
  setlist <- c(args[1])
  kvalue <- c(args[2])
  phagelist <- c(args[3])
  
  #linking phage genome accessions to their clusters
  phage <- read.table(file=phagelist, sep="|")
  colnames(phage)[4] <- c("cl")
  colnames(phage)[5] <- c("num_gm")
 
  #defining the function for maximum likelihood parameter optimization
  NLL <- function(pars, data) {
    # Extract parameters from the vector
    mu = pars[1]
    sigma = pars[2]
    # Calculate Negative Log-LIkelihood
    -sum(dnorm(x = data, mean = mu, sd = sigma, log = TRUE))
  }
  
  #summing up by samples
  for(j in 1:length(setlist)){
    
    metaset <- setlist[j]
    indir <- c(args[4])
    outdir <- c(args[5])
    
    print(paste("start to work on", metaset, "which include" , kvalue , "top hits", sep=" "))
    
    
    file <- paste(metaset,"filtered.cov.txt", sep=".")
    path <- paste(indir, kvalue,metaset,file, sep="/")
    
    print(paste("start to load hallmark abundance file", path, sep =" "))
    tab94 <- read.table(file=path)
    tab94 <- subset(tab94, tab94$V1 != "genome")

    for(i in 1:length(tab94$V1)){
      tab94$vccl[i] <- unlist(strsplit(as.character(tab94$V1[i]), "[|]"))[4]
      tab94$hallmark[i] <- unlist(strsplit(as.character(tab94$V1[i]), "[|]"))[1]
      tab94$pc[i] <- unlist(strsplit(as.character(tab94$V1[i]), "[|]"))[3]
    }
    
    uniqtab <- unique.data.frame(tab94[,6:8])
    subtmp <- subset(tab94,tab94$V2 == 0)
    testhm94 <- data.frame()
    cllist <- as.character(unique(subtmp$vccl))
    
    for(m in 1:length(cllist)) {
      
      genomes <- subset(uniqtab, uniqtab$vccl == cllist[m])$hallmark
      pc <- subset(uniqtab, uniqtab$vccl == cllist[m])$pc
      
      cov <- c()
      abd <- c()
      for(i in 1:length(genomes)){
        if( length(subset(subtmp, subtmp$hallmark == genomes[i])$V5) == 0) {
          cov[i] <- 1 
          
        } else {
          cov[i] <- 1 - subset(subtmp, subtmp$hallmark == genomes[i])$V5
        }
        
        abd[i] <- sum(subset(tab94, tab94$hallmark == genomes[i])$V2 * subset(tab94, tab94$hallmark == genomes[i])$V5)
      }
      
      tmpcl <- data.frame(genomes=genomes, pc=pc, coverage=cov, abundance=abd)
      tmpcl$cluster <- cllist[m]
      
      testhm94 <- rbind(testhm94,tmpcl)
    }
    
    testhm94 <- subset(testhm94, testhm94$coverage >= 0.25)
    
    cllist <- as.character(unique(testhm94$cluster))
    testhmcom94 <- data.frame(cluster=cllist)
    for(m in 1:length(cllist)) {
      testhmcom94$sum_abd_hm[m] <- sum(subset(testhm94, testhm94$cluster == testhmcom94$cluster[m])$abundance)
      testhmcom94$num_hallmark[m] <- length(subset(testhm94, testhm94$cluster == testhmcom94$cluster[m])$genomes)
      testhmcom94$sum_abd_hm_adj[m] <- testhmcom94$sum_abd_hm[m] / testhmcom94$num_hallmark[m]
      datat <- subset(testhm94, testhm94$cluster == testhmcom94$cluster[m])$abundance 
      
      if(sum(datat) == 0) {
        testhmcom94$nll_abd_hm[m] <- 0
        testhmcom94$nll_sd_abd_hm[m] <- 0
      }
      else if(length(datat) == 1 ) {
        testhmcom94$nll_abd_hm[m] <- datat
        testhmcom94$nll_sd_abd_hm[m] <- NA
      }
      else {
        mle = optim(par = c(mu = 10, sigma = 1), fn = NLL, data = datat)
        testhmcom94$nll_abd_hm[m] <- as.numeric(mle$par[1])
        testhmcom94$nll_sd_abd_hm[m] <- as.numeric(mle$par[2])
      }
      
    }
    
    
    testhmcom94$sum_abd_hm_adj[is.na(testhmcom94$sum_abd_hm_adj)] <- 0
    
    testhmcom94$prop_nll_abd_hm <- testhmcom94$nll_abd_hm / sum(testhmcom94$nll_abd_hm)
    testhmcom94$prop_nll_sd_abd_hm <- testhmcom94$nll_sd_abd_hm / sum(testhmcom94$nll_abd_hm)
    testhmcom94$prop_sum_abd_hm_adj <- testhmcom94$sum_abd_hm_adj / sum(testhmcom94$sum_abd_hm_adj)
    testhmcom94$log_prop_sum_abd_hm_adj <- -log10(testhmcom94$sum_abd_hm_adj / sum(testhmcom94$sum_abd_hm_adj))
    
    write.table(testhmcom94, file=paste(outdir, paste(metaset,kvalue,"txt", sep="."), sep="/"), sep="\t", quote=FALSE)
    write.table(testhm94, file=paste(outdir, paste(metaset,kvalue,"prot","txt", sep="."), sep="/"), sep="\t", quote=FALSE)
    print(paste("done with", metaset, "with parameter set" , kvalue, sep=" "))
    
  }
}
