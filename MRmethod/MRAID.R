library(MRAID)
library(parallel)


sim_cause = function(i){
    library(MRAID)
    library(dplyr)
    
    set.seed(100)
    
    # data load
    load(paste0(pathFile,"/",i,".Rdata"))
    
    X_clump = dat[(dat$p_value < 5e-8) & (dat$ld_prune == T),]

    
    if(dim(X_clump)[1] > 3){
        #LD_matrix = 
        #    ieugwasr::ld_matrix(variants = X_clump$rsid,
        #                        plink_bin = "/home/tb/tools/GWAS//plink",
        #                        bfile="/home/tb/ref/1kg.v3//EUR")
        
        # MRAID fit
        res = MRAID(Zscore_1 = X_clump$beta_hat_1 / X_clump$seb1, Zscore_2 = X_clump$beta_hat_2/X_clump$seb2,
                   Sigma1sin = diag(rep(1,dim(X_clump)[1])),Sigma2sin = diag(rep(1,dim(X_clump)[1])),
                   samplen1 = as.numeric(n1),samplen2 = as.numeric(n2))
                   
        # save results
        P = res$causal_pvalue
        eff = res$causal_effect
        correlated_eff = res$correlated_pleiotropy_effect
        
        cat(paste(i,dim(X_clump)[1],P,eff,correlated_eff,"\n",sep="\t"), file=saveFile,append=TRUE)
    }else{
    
        # save results
        P = 0
        eff = 0
        
        cat(paste(i,dim(X_clump)[1],P,eff,correlated_eff,"\n",sep="\t"), file=saveFile,append=TRUE)
        
        # print(paste0(folder,i,"completed"))
    }

    

}


simulationTimes = 100
Paras =  readLines("./MRmethod/allParas.txt")
path = Paras[1]
savePath = Paras[2]

folders = readLines("./MRmethod/all_runFiles.txt")
  
  

for(folder in folders){
  
  if(!file.exists(paste(savePath,folder,sep="/"))){
    dir.create(paste(savePath,folder,sep="/"))
  }
  
  pathFile = paste(path,folder,sep="/")
  saveFile = paste(savePath,folder,"MRAID.txt",sep="/")
  file.create(paste(saveFile,sep="/"))
  cat(paste("simNum","IVsNum","P","eff","correlated_eff\n",sep="\t"), file=saveFile)
  
  print(pathFile)
  print(saveFile)
  
  n1 = strsplit(folder,"_")[[1]][7]
  n2 = strsplit(folder,"_")[[1]][8]
  
  cl <- makeCluster(10)
  clusterExport(cl,"saveFile",envir = environment())
  clusterExport(cl,"pathFile",envir = environment())
  clusterExport(cl,"folder",envir = environment())
  clusterExport(cl,"n1",envir = environment())
  clusterExport(cl,"n2",envir = environment())

  parLapply(cl,seq(simulationTimes),sim_cause)

  stopCluster(cl) # ?رռ?Ⱥ
}