library(parallel) # ???߳?


sim_cause = function(i){
    library(stringr)
    library(PCMR)
    library(weights)
    
    # data load
    load(paste0(pathFile,"/",i,".Rdata"))
    
    set.seed(0)    
            
    # PCMR process-------------------------- 
    # evaluate uncorrelated pleiotropy
    X_clump1 = dat[sample(seq(dim(dat)[1]),10000),]
  
    init = PCMR_initEst(X_clump1$beta_hat_1,X_clump1$seb1,
                          X_clump1$beta_hat_2,X_clump1$seb2)
    
    # IVs 
    X_clump = dat[(dat$p_value < 5e-8) & (dat$ld_prune == T),]  
    
    if(dim(X_clump)[1] > 3){
        # PCMR fit
        num_gamma = 2
        # model 1 : gSigma2 T ; sigma2 T
        result_1 =  PCMR(X_clump$beta_hat_1, X_clump$seb1,
                         X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="1",
                         # sigma2 = init$sigma2,rho=init$rho,isIntact=T)
                         sigma2 = init$sigma2,isIntact=F)
                                                  
        result_1 = PCMR_cEst(result_1,ref_beta_outcome = X_clump1$beta_hat_2,ref_se_outcome = X_clump1$seb2,cores=12) # estimate the factor c
        result_1 = PCMR_testCausal_bootstrap(result_1,cores=12) # bootstrapping to estimate D_HVP
        result_1 = PCMR_testCorPlei(result_1) # calculate Pvalue of heterogeneity
    
        result_1 = PCMR_testCausal(result_1)          
        
        # save results
        P_1 = result_1$Pvalue
        eff_1 = result_1$effect
        
        eff_mean_minClass = result_1$bootstrap$mean_minClass 
        sd_minClass = result_1$bootstrap$sd_minClass 
        P_minClass_mean  = result_1$bootstrap$P_minClass_mean 
        
        eff_mean_maxClass = result_1$bootstrap$mean_maxClass
        sd_maxClass = result_1$bootstrap$sd_maxClass
        P_maxClass_mean  = result_1$bootstrap$P_maxClass_mean 
        
        correct_factor = result_random$c
        CHVP_test_correct = result_1$CHVP_test
        
        AIC_1 = result_1$AIC
        gamma_1 = str_c(result_1$gamma,collapse = "_")
        prop_1 = str_c(result_1$pi_gamma,collapse = "_")
        
        cat(paste(i,dim(X_clump)[1],
                  P_random,eff_random,AIC_random,gamma_random,prop_random,
                  eff_mean_minClass,sd_minClass,P_minClass_mean,
                  eff_mean_maxClass,sd_maxClass,P_maxClass_mean,
                  correct_factor,CHVP_test_correct,"\n",sep="\t"), file=saveFile,append=TRUE)
        
        # print(paste0(folder,i,"completed"))
    }else{
        cat(paste(i,dim(X_clump)[1],rep(0,13),"\n",sep="\t"), file=saveFile,append=TRUE)
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
  saveFile = paste(savePath,folder,"PCMR(n=2)_intact_bootstrap_interval.txt",sep="/")
  
  SimulationNums = seq(simulationTimes)
  if(!file.exists(saveFile)){
    file.create(paste(saveFile,sep="/"))
    cat(paste("simNum","IVsNum",
              "P_random","eff_random","AIC_random","gamma_random","prop_random",
              "eff_mean_minClass","sd_minClass","P_minClass_mean",
              "eff_mean_maxClass","sd_maxClass","P_maxClass_mean",
              "correct_factor","CHVP_test_correct\n",sep="\t"), file=saveFile)
  }else{
    saveTable = read.table(saveFile,header=1)
    SimulationNums = setdiff(SimulationNums,saveTable$simNum)
  }
  
  
  print(pathFile)
  print(saveFile)
  
  cl <- makeCluster(10)
  clusterExport(cl,"saveFile",envir = environment())
  clusterExport(cl,"pathFile",envir = environment())
  clusterExport(cl,"folder",envir = environment())
  
  clusterApplyLB(cl,SimulationNums,sim_cause)

  stopCluster(cl) 
}