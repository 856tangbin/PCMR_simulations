library(cause)
library(parallel) # ???߳?

saveSvg = function(result,S,model){
    svg(
        filename = paste0("./Application/Psychiatric_disorders_analysis/results/PCMR(n=1)_intact/",S[1],"_",S[2],"model",model,".svg"), # ?ļ?????
        width = 7,           # ??
        height = 7,          # ??
    )              # ?ֱ???
    PCMR_plot(result)
    dev.off()
}

sim_cause = function(i){
    library(cause)
    library(readr)
    library(dplyr)
    library(PCMR)
    
    Dir = "./Application/Psychiatric_disorders_analysis/trait/"
    trait = read.csv(paste0(Dir,"trait.csv"))
    tests = read.csv(paste0(Dir,"common disease.csv")) 
    
    # data load
    test = tests$Traits[i]
    S =  stringr::str_split(test," -> ")[[1]]
    risk = trait[trait$Abbreviation == tolower(S[1]),]
    disease = trait[trait$Abbreviation == tolower(S[2]),]
    
    # load data
    risk_data = read_delim(paste0(Dir,risk$Filename),delim=list(space=" ",tab="\t")[[risk$delimeter]],skip=risk$skip)
    disease_data = read_delim(paste0(Dir,disease$Filename),delim=list(space=" ",tab="\t")[[disease$delimeter]],skip=disease$skip)
    
    # combine gwas summary by cause
    risk_cols = risk[,c("rsid","beta","se","A1","A2","P.value")]
    disease_cols = disease[,c("rsid","beta","se","A1","A2","P.value")]
    
    # numeric
    for(ind in c("beta","se","P.value")){
        risk_data[[risk_cols[[ind]]]] = as.numeric(risk_data[[risk_cols[[ind]]]])
        disease_data[[disease_cols[[ind]]]] = as.numeric(disease_data[[disease_cols[[ind]]]])
    }
    
    X <- gwas_merge(risk_data, disease_data, snp_name_cols = c(risk_cols[[1]], disease_cols[[1]]), 
                    beta_hat_cols = c(risk_cols[[2]], disease_cols[[2]]), 
                    se_cols = c(risk_cols[[3]], disease_cols[[3]]), 
                    A1_cols = c(risk_cols[[4]], disease_cols[[4]]), 
                    A2_cols = c(risk_cols[[5]], disease_cols[[5]]), 
                    pval_cols = c(risk_cols[[6]], disease_cols[[6]]))
    
    r2_thresh = 0.1
    pval_thresh = 5e-8
    
    # LD pruning
    X_clump <- X[X$p1 < pval_thresh,] %>%
        rename(rsid = snp,
               pval = p1) %>%
        ieugwasr::ld_clump(dat = .,
                           clump_r2 = r2_thresh,
                           clump_p = pval_thresh,
                           plink_bin = "/home/tb/tools/GWAS//plink",
                           bfile="/home/tb/ref/1kg.v3//EUR")
    
    X_ = X[X$p1 > 0.5,]
    X_clump1 <- X_[sample(seq(dim(X_)[1]),100000),] %>%
        rename(rsid = snp,
               pval = p1) %>%
        ieugwasr::ld_clump(dat = .,
                           clump_r2 = r2_thresh,
                           clump_p = 1,
                           plink_bin = "/home/tb/tools/GWAS//plink",
                           bfile="/home/tb/ref/1kg.v3/EUR")
    
    # init calculate
    init = PCMR_initEst(X_clump1$beta_hat_1,X_clump1$seb1,
                        X_clump1$beta_hat_2,X_clump1$seb2)

    # model fitting
    num_gamma = 1
    
    # model 1 : gSigma2 T ; sigma2 T
    result_1 = PCMR(X_clump$beta_hat_1, X_clump$seb1,
                    X_clump$beta_hat_2,X_clump$seb2,num_gamma = num_gamma,model="1",
                    sigma2 = init$sigma2,rho=init$rho,isIntact = T)  
    saveSvg(result_1,S,1)

    # model 2 : gSigma2 F ; sigma2 T
    result_2 =  PCMR(X_clump$beta_hat_1, X_clump$seb1,
                    X_clump$beta_hat_2,X_clump$seb2,num_gamma = num_gamma,model="2",
                    sigma2 = init$sigma2,rho=init$rho,isIntact = T)
    saveSvg(result_2,S,2)
    
    # model 3 : gSigma2 T ; sigma2 F
    result_3 =  PCMR(X_clump$beta_hat_1, X_clump$seb1,
                    X_clump$beta_hat_2,X_clump$seb2,num_gamma = num_gamma,model="3",
                    sigma2 = init$sigma2,rho=init$rho,isIntact = T)  
    saveSvg(result_3,S,3)
    
    # model 4 : gSigma2 F ; sigma2 F
    result_4 =  PCMR(X_clump$beta_hat_1, X_clump$seb1,
                    X_clump$beta_hat_2,X_clump$seb2,num_gamma = num_gamma,model="4",
                    sigma2 = init$sigma2,rho=init$rho,isIntact = T)
    saveSvg(result_4,S,4)
    
    # save results
    P_1 = result_1$Pvalue
    eff_1 = result_1$effect
    AIC_1 = result_1$AIC
  
    P_2 = result_2$Pvalue
    eff_2 = result_2$effect
    AIC_2 = result_2$AIC
        
    P_3 = result_3$Pvalue
    eff_3 = result_3$effect
    AIC_3 = result_3$AIC
    
    P_4 = result_4$Pvalue
    eff_4 = result_4$effect
    AIC_4 = result_4$AIC
    
    cat(paste(risk$Abbreviation,disease$Abbreviation,dim(X_clump)[1],P_1,eff_1,AIC_1,P_2,eff_2,AIC_2,P_3,eff_3,AIC_3,P_4,eff_4,AIC_4,"\n",sep="\t"), file=saveFile,append=TRUE)
}


Dir = "./Application/Psychiatric_disorders_analysis/trait/"
savePath = "./Application/Psychiatric_disorders_analysis/results/"
tests = read.csv(paste0(Dir,"common disease.csv")) 

saveFile = paste(savePath,"PCMR(n=1)_intact.txt",sep="/")
file.create(paste(saveFile,sep="/"))
cat(paste("risk","disease","IVs","P_1","eff_1","AIC_1","P_2","eff_2","AIC_2","P_3","eff_3","AIC_3","P_4","eff_4","AIC_4\n",sep="\t"), file=saveFile)

cl <- makeCluster(48)
clusterExport(cl,"saveFile",envir = environment())
clusterExport(cl,"saveSvg",envir = environment())
parLapply(cl,seq(dim(tests)[1]),sim_cause)

stopCluster(cl) # ?رռ?Ⱥ
