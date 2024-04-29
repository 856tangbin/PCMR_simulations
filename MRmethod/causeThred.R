library(cause)
library(parallel) 


sim_cause = function(i){
    library(cause)
    library(dplyr)
    
    set.seed(100)
    
    # data load
    load(paste0(pathFile,"/",i,".Rdata"))
    
    # CAUSE process--------------------------
    # Calculate nuisance parameters
    varlist <- with(dat, sample(snp, size=100000, replace=FALSE))
    params <- est_cause_params(dat, varlist)
    
    X_clump = dat[(dat$p_value < 1e-3) & (dat$ld_prune == T),]

          
    if(dim(X_clump)[1] > 5){
        # CAUSE fit
        res <- cause(X=dat, variants = X_clump$snp, param_ests = params) 
        result = summary(res, ci_size=0.95)
        
        # save results
        eff = result$quants[[2]][1,1]
        q_median_causal = result$quants[[2]][1,3]
        eta_median_causal = result$quants[[2]][1,2]
        
        q_median_sharing = result$quants[[1]][1,3]
        eta_median_sharing = result$quants[[1]][1,2]
        P = result$p
        
        rho = params$rho
        
        
        # P = result[3]
        # eff = summary(res, ci_size=0.95)[1]$quants[[2]][1,1]
        
         cat(paste(i,dim(X_clump)[1], eff,	q_median_causal,	eta_median_causal,	q_median_sharing,	eta_median_sharing,	P,	rho,"\n",sep="\t"), file=saveFile,append=TRUE)
    }else{
    
        
        # save results
        eff = 0
        q_median_causal =0
        eta_median_causal = 0
        
        q_median_sharing = 0
        eta_median_sharing = 0
        P = 0
        
        cat(paste(i,dim(X_clump)[1], eff,	q_median_causal,	eta_median_causal,	q_median_sharing,	eta_median_sharing,	P,	rho,"\n",sep="\t"), file=saveFile,append=TRUE)
        
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
  saveFile = paste(savePath,folder,"causeThred.txt",sep="/")
  file.create(paste(saveFile,sep="/"))
  cat(paste("simNum","IVsNum","eff","q_median_causal","eta_median_causal","q_median_sharing","eta_median_sharing","P","rho\n",sep="\t"), file=saveFile)
  
  print(pathFile)
  print(saveFile)
  
  cl <- makeCluster(8)
  clusterExport(cl,"saveFile",envir = environment())
  clusterExport(cl,"pathFile",envir = environment())
  clusterExport(cl,"folder",envir = environment())
  
  parLapply(cl,seq(simulationTimes),sim_cause)

  stopCluster(cl) 
}