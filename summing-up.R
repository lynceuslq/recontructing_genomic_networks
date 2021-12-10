#!/bin/Rscript

setlist <- read.table("/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/testongenomes/abd_dist_model_co_split/MYSAMPLIST")$V1
phagelist <- c("/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/rm_overlap/gpdvcmorethan3.clustered.len.overlappingremoved.phage.list")

#linking phage genome accessions to their clusters
phage <- read.table(file=phagelist, sep="|")
for(i in 1:length(phage$V1)){
  phage$cl[i] <-  unlist(strsplit(as.character(phage$V4[i]), ":"))[1]
  phage$num_gm[i] <- unlist(strsplit(as.character(phage$V4[i]), ":"))[2]
}

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
  indir <- c("/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/testongenomes/abd_dist_model_co_split")
  outdir <- c("/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/testongenomes/abd_dist_model_co_split/by_gm_cov10")
  simdir <- c("/jdfssz1/ST_HEALTH/P20Z10200N0206/liqian6/meta-sim/abd_dist_model_testsets/abd_pack")
  
  print(paste("start to work on", metaset, sep=" "))
  
  file <- paste(metaset,"cov.txt", sep="_")
  path <- paste(indir,metaset,file, sep="/")
  
  print(paste("start to load genome avrg depth file", path, sep =" "))
  testbygm <- read.table(file=path)
  testbygm <- subset(testbygm, testbygm$V1 != "genome")
  
  file <- paste(metaset, "abundance.txt", sep="_")
  path <- paste(simdir, file, sep="/")
  print(paste("start to load simulated abundance file", path, sep =" "))
  
  abd_sum <- read.table(file=path)
  abd_sum$cl <- phage$cl[match(abd_sum$V1, phage$V1)]
  
  cllist <- as.character(unique(phage$cl))
  subtmp <- subset(testbygm,testbygm$V2 == 0)
  subtmp <- subset(subtmp,subtmp$V5 <= 0.9)
  test94 <- data.frame()  

  for(m in 1:length(cllist)) {
    genomes <- as.character(subset(phage, phage$cl == cllist[m])$V1)
    
    cov <- c()
    abd <- c()
    for(i in 1:length(genomes)){
      cov[i] <- 1 - subtmp[match(genomes[i], subtmp$V1), ]$V5
    }
    
    for(i in 1:length(genomes)){
      abd[i] <- sum(subset(testbygm, testbygm$V1 == genomes[i])$V2 * subset(testbygm, testbygm$V1 == genomes[i])$V5)
    }
    
    tmpcl <- data.frame(genomes=genomes, coverage=cov, abundance=abd)
    tmpcl$cluster <- cllist[m]
    
    test94 <- rbind(test94,tmpcl)
  }
  
  write.table(test94, file=paste(outdir, paste(metaset,"genomes","txt", sep="."), sep="/"), sep="\t", quote=FALSE)
  
  cllist <- as.character(unique(test94$cluster))
  testhmcom94 <- data.frame(cluster=cllist)
  for(m in 1:length(cllist)) {
    # testhmcom94$sum_abd_hm[m] <- sum(subset(testhm94, testhm94$cluster == testhmcom94$cluster[m])$abundance)
    # testhmcom94$num_hallmark[m] <- length(subset(testhm94, testhm94$cluster == testhmcom94$cluster[m])$genomes)
    # testhmcom94$sum_abd_hm_adj[m] <- testhmcom94$sum_abd_hm[m] / testhmcom94$num_hallmark[m]
    testhmcom94$ref_sum[m] <- sum(subset(abd_sum, abd_sum$cl == testhmcom94$cluster[m])$V2)
    
    testhmcom94$sum_abd_gm[m] <- sum(subset(test94, test94$cluster == testhmcom94$cluster[m])$abundance) 
    testhmcom94$num_genome[m] <- length(subset(test94, test94$cluster == testhmcom94$cluster[m])$genomes )
    testhmcom94$sum_abd_gm_adj[m] <- testhmcom94$sum_abd_gm[m] / testhmcom94$num_genome[m]
    testhmcom94$ref_sum[m] <- sum(subset(abd_sum, abd_sum$cl == testhmcom94$cluster[m])$V2)
#    datat <- subset(test94, test94$cluster == testhmcom94$cluster[m])$abundance    
    
#    if(sum(datat) == 0) {
#      testhmcom94$nll_abd_gm[m] <- 0
#      testhmcom94$nll_sd_abd_gm[m] <- 0
#    }
#    else if(length(datat) == 1 ) {
#      testhmcom94$nll_abd_gm[m] <- datat
#      testhmcom94$nll_sd_abd_gm[m] <- NA
#    }
#    else {
#      mle = optim(par = c(mu = 10, sigma = 1), fn = NLL, data = datat)
#      testhmcom94$nll_abd_gm[m] <- as.numeric(mle$par[1])
#      testhmcom94$nll_sd_abd_gm[m] <- as.numeric(mle$par[2])
#    }
    
  }
  
  
  testhmcom94$sum_abd_gm_adj[is.na(testhmcom94$sum_abd_gm_adj)] <- 0
  
  testhmcom94$adj_ref_sum <- testhmcom94$ref_sum / sum(testhmcom94$ref_sum)
  testhmcom94$prop_sum_abd_gm  <- testhmcom94$sum_abd_gm / sum(testhmcom94$sum_abd_gm)
#  testhmcom94$prop_nll_abd_gm <- testhmcom94$nll_abd_gm / sum(testhmcom94$nll_abd_gm)
#  testhmcom94$prop_nll_sd_abd_gm <- testhmcom94$nll_sd_abd_gm / sum(testhmcom94$nll_abd_gm)
  testhmcom94$prop_sum_abd_gm_adj <- testhmcom94$sum_abd_gm_adj / sum(testhmcom94$sum_abd_gm_adj)
  testhmcom94$log_prop_sum_abd_gm_adj <- -log10(testhmcom94$sum_abd_gm_adj / sum(testhmcom94$sum_abd_gm_adj))
  
  write.table(testhmcom94, file=paste(outdir, paste(metaset,"genomic_abd","txt", sep="."), sep="/"), sep="\t", quote=FALSE)
  print(paste("done with", metaset, sep=" "))
  
}

