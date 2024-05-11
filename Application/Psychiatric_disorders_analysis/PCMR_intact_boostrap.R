library(cause)
library(parallel) 

sim_cause = function(i){
    library(cause)
    library(readr)
    library(dplyr)
    library(PCMR)
    library(weights)
    
    library(stringr)
    
    set.seed(0)
    
    Dir = "./Application/Common_disease_analysis/trait/"
    setwd(Dir)
    trait = read.csv(paste0(Dir,"trait.csv"))
    tests = read.csv(paste0(Dir,"common disease.csv")) 
    
    # data load
    test = tests$Traits[i]
    S =  stringr::str_split(test," -> ")[[1]]
    
    print(S)
    
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
    
    X1 = X[X$p1 > 0.5,]
    X_clump1 <- X1[sample(seq(dim(X1)[1]),100000),] %>%
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
    num_gamma = 2
    
    # model 1 : gSigma2 T ; sigma2 T
    result_1 =  PCMR(X_clump$beta_hat_1, X_clump$seb1,
                     X_clump$beta_hat_2,X_clump$seb2,num_gamma = 2,model="1",
                     sigma2 = init$sigma2,rho=init$rho,isIntact=T)
    
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
    
    cat(paste(risk$Abbreviation,disease$Abbreviation,dim(X_clump)[1],
              P_random,eff_random,AIC_random,gamma_random,prop_random,
              eff_mean_minClass,sd_minClass,P_minClass_mean,
              eff_mean_maxClass,sd_maxClass,P_maxClass_mean,
              correct_factor,CHVP_test_correct,"\n",sep="\t"), file=saveFile,append=TRUE)
    print(paste0(S," end"))
}


Dir = "./Application/Psychiatric_disorders_analysis/trait/"
savePath = "./Application/Psychiatric_disorders_analysis/results/"
tests = read.csv(paste0(Dir,"common disease.csv")) 

saveFile = paste(savePath,"PCMR(n=2)_intact_boostrap_interval.txt",sep="/")
file.create(paste(saveFile,sep="/"))
cat(paste("risk","disease","IVs",
          "P_random","eff_random","AIC_random","gamma_random","prop_random",
          "eff_mean_minClass","sd_minClass","P_minClass_mean",
          "eff_mean_maxClass","sd_maxClass","P_maxClass_mean",
          "correct_factor","CHVP_test_correct\n",sep="\t"), file=saveFile)

for(i in seq(dim(tests)[1])){
    sim_cause(i)
}


