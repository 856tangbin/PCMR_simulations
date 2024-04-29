library(MRPRESSO)
library(parallel) 


sim_cause = function(i){
    library(MRPRESSO)
    library(dplyr)
    
    # data load
    load(paste0(pathFile,"/",i,".Rdata"))
    
    X_clump = dat[(dat$p_value < 5e-8) & (dat$ld_prune == T),]

      
    if(dim(X_clump)[1] > 3){
        # mrpresso fit
        res = mr_presso(BetaOutcome = "beta_hat_2", BetaExposure = "beta_hat_1", SdOutcome = "seb2", SdExposure = "seb1", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = X_clump, NbDistribution = 1000,  SignifThreshold = 0.05)
        
        
        
        # save results
        plei_test = res$`MR-PRESSO results`$`Global Test`$Pvalue
        
        if(is.character(plei_test)){
            plei_test = 0
        }
        
        P = res$`Main MR results`$`P-value`[2]
        eff = res$`Main MR results`$`Causal Estimate`[2]
        
        if(is.na(P)){
            P = res$`Main MR results`$`P-value`[1]
            eff = res$`Main MR results`$`Causal Estimate`[1]
        }
        
        cat(paste(i,dim(X_clump)[1],plei_test,P,eff,"\n",sep="\t"), file=saveFile,append=TRUE)
    }else{
    
        # save results
        P = 0
        eff = 0
        cat(paste(i,dim(X_clump)[1],plei_test,P,eff,"\n",sep="\t"), file=saveFile,append=TRUE)
    }
    

    
    
    # print(paste0(folder,i,"completed"))
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
  saveFile = paste(savePath,folder,"mrpresso.txt",sep="/")
  file.create(paste(saveFile,sep="/"))
  cat(paste("simNum","IVsNum","plei_test","P","eff\n",sep="\t"), file=saveFile)
  
  print(pathFile)
  print(saveFile)
  
  cl <- makeCluster(50)
  clusterExport(cl,"saveFile",envir = environment())
  clusterExport(cl,"pathFile",envir = environment())
  clusterExport(cl,"folder",envir = environment())
  
  parLapply(cl,seq(simulationTimes),sim_cause)

  stopCluster(cl)
}