library(stringr)
library(cause)
library(causeSims)
library(furrr)
library(RcppZiggurat)
library(MendelianRandomization)


Path = "../"
savePath = paste0(Path,"causeSimData/")

snps <- readRDS(paste0(Path,"data/chr19_snpdata_hm3only.RDS")) # download through https://github.com/jean997/causeSims
evd_list <- readRDS(paste0(Path,"data/evd_list_chr19_hm3.RDS"))

simData_ex = function(gamma,q,eta,simulationTimes = 100,n1=40000,n2=40000,neffect1= 1000,neffect2= 1000,h1= 0.25,h2= 0.25,seed = 1){
  
  tau = gamma**2
  omega = q*eta**2 

  if(gamma == 0){
    folder = paste("H0",gamma,tau,q,eta,omega,n1,n2,neffect1,neffect2,h1,h2,sep="_")  
  }else{
    folder = paste("H1",gamma,tau,q,eta,omega,n1,n2,neffect1,neffect2,h1,h2,sep="_")  
  }
  
  
  if(!file.exists(paste(savePath,folder,sep="/"))){
    dir.create(paste(savePath,folder,sep="/"))
  }

  for(. in seq(simulationTimes)){
    
    filename = paste0(savePath,folder,"/",. ,".Rdata")
    if(file.exists(filename) & (file.size(filename) > 1024 * 1000)){
      next()
    }
    dat <- sum_stats(snps, evd_list,
                     n_copies = 30,
                     n1 = n1, n2=n2,
                     h1=h1, h2=h2,
                     neffect1 = neffect1, neffect2 = neffect2,
                     gamma = gamma, eta = eta, q = q,
                     cores = 1, ld_prune_pval_thresh = 0.01,
                     r2_thresh = 0.1,seed=as.integer(runif(1,1,1e8)))

    save(dat,file = filename)
  }
}

