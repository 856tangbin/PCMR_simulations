library(cause)
library(parallel) 


sim_cause = function(i){
    library(MendelianRandomization)
    library(dplyr)
    
    # data load
    load(paste0(pathFile,"/",i,".Rdata"))
    
    X_clump = dat[(dat$p_value < 5e-8) & (dat$ld_prune == T),]
    
    if(dim(X_clump)[1] > 3){
        # MendelianRandomization -------------------------
        mr_data = mr_input(bx = X_clump$beta_hat_1, bxse = X_clump$seb1,
                           by = X_clump$beta_hat_2, byse = X_clump$seb2)
        
        result = mr_mbe(mr_data)
        P_mbe = result$Pvalue
        eff_mbe = result$Estimate
        
        result = mr_egger(mr_data)
        P_egger = result$Pvalue.Est
        eff_egger = result$Estimate
        
        result = mr_ivw(mr_data)
        P_ivw = result$Pvalue
        eff_ivw = result$Estimate
        
        result = mr_median(mr_data)
        P_weightedMedian = result$Pvalue
        eff_weightedMedian = result$Estimate
        
        # save results
        cat(paste(i,dim(X_clump)[1],P_mbe,eff_mbe,P_egger,eff_egger,P_ivw,eff_ivw,P_weightedMedian,eff_weightedMedian,"\n",sep="\t"), file=saveFile,append=TRUE)
    }else{
        # save results
        cat(paste(i,dim(X_clump)[1],0,0,0,0,0,0,0,0,"\n",sep="\t"), file=saveFile,append=TRUE)
        
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
  saveFile = paste(savePath,folder,"mendelianRandomization.txt",sep="/")
  file.create(paste(saveFile,sep="/"))
  cat(paste("simNum","IVsNum","P_mbe","eff_mbe","P_egger","eff_egger","P_ivw","eff_ivw","P_weightedMedian","eff_weightedMedian\n",sep="\t"), file=saveFile)
  
  print(pathFile)
  print(saveFile)
  
  cl <- makeCluster(20)
  clusterExport(cl,"saveFile",envir = environment())
  clusterExport(cl,"pathFile",envir = environment())
  clusterExport(cl,"folder",envir = environment())
  
  parLapply(cl,seq(simulationTimes),sim_cause)

  stopCluster(cl) 
}