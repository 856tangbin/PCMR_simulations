library(parallel) 


sim_cause = function(i){
    library(stringr)
    library(PCMR)
    library(weights)
    
    # data load
    load(paste0(pathFile,"/",i,".Rdata"))
    
    # PCMR process--------------------------
    # evaluate uncorrelated pleiotropy
    X_clump1 = dat[sample(seq(dim(dat)[1]),10000),]
  
    init = PCMR_initEst(X_clump1$beta_hat_1, X_clump1$seb1,
                        X_clump1$beta_hat_2,X_clump1$seb2)
    
    # IVs 
    X_clump = dat[(dat$p_value < 5e-8) & (dat$ld_prune == T),]  
    
    if(dim(X_clump)[1] > 3){
        # PCMR fit
        num_gamma = 2
        # model 1 : gSigma2 T ; sigma2 T
        result_1 =  PCMR(X_clump$beta_hat_1, X_clump$seb1,
                         X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="1",
                         sigma2 = init$sigma2,rho=init$rho,isIntact=T)
    
        # model 2 : gSigma2 F ; sigma2 T
        result_2 =  PCMR(X_clump$beta_hat_1, X_clump$seb1,
                         X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="2",
                         sigma2 = init$sigma2,rho=init$rho,isIntact=T)
       
        # save results
        P_1 = result_1$Pvalue
        eff_1 = result_1$effect
        AIC_1 = result_1$AIC
        gamma_1 = str_c(result_1$gamma,collapse = "_")
        prop_1 = str_c(result_1$pi_gamma,collapse = "_")
      
        P_2 = result_2$Pvalue
        eff_2 = result_2$effect
        AIC_2 = result_2$AIC
        gamma_2 = str_c(result_2$gamma,collapse = "_")
        prop_2 = str_c(result_2$pi_gamma,collapse = "_")
        
        cat(paste(i,dim(X_clump)[1],
                  P_random,eff_random,AIC_random,gamma_random,prop_random,
                  P_fixed,eff_fixed,AIC_fixed,gamma_fixed,prop_fixed,"\n",sep="\t"), file=saveFile,append=TRUE)
        
        # print(paste0(folder,i,"completed"))
    }else{
        cat(paste(i,dim(X_clump)[1],
                  0,0,0,0,0,0,
                  0,0,0,0,0,0,"\n",sep="\t"), file=saveFile,append=TRUE)
    
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
  saveFile = paste(savePath,folder,"PCMR(n=2)_intact.txt",sep="/")
  file.create(paste(saveFile,sep="/"))
  cat(paste("simNum","IVsNum",
            "P_random","eff_random","AIC_random","gamma_random","prop_random",
            "P_fixed","eff_fixed","AIC_fixed","gamma_fixed","prop_fixed\n",sep="\t"), file=saveFile)
  
  print(pathFile)
  print(saveFile)
  
  cl <- makeCluster(80)
  clusterExport(cl,"saveFile",envir = environment())
  clusterExport(cl,"pathFile",envir = environment())
  clusterExport(cl,"folder",envir = environment())
  
  parLapply(cl,seq(simulationTimes),sim_cause)

  stopCluster(cl) 
}