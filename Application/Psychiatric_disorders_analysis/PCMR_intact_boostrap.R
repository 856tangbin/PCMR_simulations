library(cause)
library(parallel) # ???߳?


boot = function(.,result){
    library(PCMR)
    
    size = length(result$Paras$beta_ex)
    Ind = sample(seq(size),size = size,replace = T)
    
    Paras = result$Paras
    Paras$beta_ex = result$Paras$beta_ex[Ind]
    Paras$beta_ex_se = result$Paras$beta_ex_se[Ind]
    Paras$beta_out = result$Paras$beta_out[Ind]
    Paras$beta_out_se = result$Paras$beta_out_se[Ind]
    results0 = pleiClassify(Paras)
    
    return(sort(results0$gamma))
}

boot_sample = function(.,result,beta_outcome,se_outcome){
    library(PCMR)
    
    result_sample = result
    
    result_sample$Paras$beta_out = beta_outcome[Seqs[[.]]] + result$est1 * rnorm(length(result$Paras$beta_ex),result$Paras$beta_ex,result$Paras$beta_ex_se)
    # result_sample$Paras$beta_out_se = sqrt(se_outcome[Seqs[[.]]]**2 + (result$est1 * result$Paras$beta_ex_se)**2)
    result_sample$Paras$beta_out_se = se_outcome[Seqs[[.]]]
    
    return(boot(.,result_sample))
}

PCMR_correct = function(result,beta_outcome,se_outcome,samples=100,sample_boot = 30,cores=80){


    result$est1 = PCMR(result$Paras$beta_ex, result$Paras$beta_ex_se,
                         result$Paras$beta_out,result$Paras$beta_out_se,num_gamma = 1,model="1")$gamma
    
    
    Seqs = list()
    for(i in seq(samples*sample_boot)){
        
        if(i %% sample_boot == 1){
            Seqs[[i]] = sample(seq(length(beta_outcome)),(length(result$Paras$beta_ex)))
        }else{
            Seqs[[i]] = Seqs[[i-1]]
        }
    }
    
    sub_cl <- makeCluster(cores)

    clusterExport(sub_cl,"Seqs",envir = environment())  
    clusterExport(sub_cl,"boot",envir = environment())  
    Cup = clusterApplyLB(sub_cl,seq(samples*sample_boot),boot_sample,result,beta_outcome,se_outcome)
    Cup_gamma = matrix(unlist(Cup),ncol=2,byrow=T)

    stopCluster(sub_cl) # ?رռ?Ⱥ
    
    Cup_chi2 = c()
    for(i in seq(samples)){
        
        Cup_gamma_sub = Cup_gamma[seq((i-1)*sample_boot+1,i*sample_boot),]
        Cup_chi2 = c(Cup_chi2,(mean(Cup_gamma_sub[,1]) - mean(Cup_gamma_sub[,2]))^2/( var(Cup_gamma_sub[,1]) +  var(Cup_gamma_sub[,2])))
    }
    
    return(Cup_chi2)
}

PCMR_testCausal_bootstrap = function(result,samples=1000,cores=80){
    
    sub_cl <- makeCluster(cores)
    clusterExport(sub_cl,"result",envir = environment())
    
    Cup = parSapply(sub_cl,seq(samples),boot,result,simplify = TRUE)
    stopCluster(sub_cl) # ?رռ?Ⱥ
    
    Cup_gamma = matrix(unlist(Cup),ncol=2,byrow=T)
    
    result$mean_minClass = mean(Cup_gamma[,1])
    result$mean_maxClass = mean(Cup_gamma[,2])
    
    result$sd_minClass = sd(Cup_gamma[,1])
    result$sd_maxClass = sd(Cup_gamma[,2])
    
    result$min_interval = stringr::str_c(quantile(Cup_gamma[,1],c(0.025,0.975)),collapse = "_")
    result$max_interval = stringr::str_c(quantile(Cup_gamma[,2],c(0.025,0.975)),collapse = "_")
    
    result$CHVP_test = pchisq((result$mean_minClass - result$mean_maxClass)^2/( result$sd_minClass^2 +  result$sd_maxClass^2),1,lower.tail = F)
    
    result$P_minClass = pt(abs(min(result$gamma)/result$sd_minClass), df=length(result$Paras$beta_ex)-1, lower.tail=F)*2
    result$P_maxClass = pt(abs(max(result$gamma)/result$sd_maxClass), df=length(result$Paras$beta_ex)-1, lower.tail=F)*2
    
    result$P_minClass_mean = pt(abs(mean(Cup_gamma[,1])/result$sd_minClass), df=length(result$Paras$beta_ex)-1, lower.tail=F)*2
    result$P_maxClass_mean = pt(abs(mean(Cup_gamma[,2])/result$sd_maxClass), df=length(result$Paras$beta_ex)-1, lower.tail=F)*2
    
    result$P_minClass_permut = 1 - abs(sum(Cup_gamma[,1] > 0) - sum(Cup_gamma[,1] < 0))/samples
    result$P_maxClass_permut = 1 - abs(sum(Cup_gamma[,2] > 0) - sum(Cup_gamma[,2] < 0))/samples
    
    return(result)
}

sim_cause = function(i){
    library(cause)
    library(readr)
    library(dplyr)
    library(PCMR)
    library(weights)
    
    library(stringr)
    
    set.seed(0)
    
    Dir = "./Application/Psychiatric_disorders_analysis/trait/"
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
    result_1 = PCMR(X_clump$beta_hat_1, X_clump$seb1,
                    X_clump$beta_hat_2,X_clump$seb2,num_gamma = num_gamma,model="1",
                    sigma2 = init$sigma2,rho=init$rho,isIntact = T)  
    result_1 = PCMR_testCorPlei(result_1)
    
    result_1 = PCMR_testCausal_bootstrap(result_1)
    
    # save results
    P_1 = result_1$Pvalue
    eff_1 = result_1$effect
    
    
    mean_minClass = result_1$mean_minClass
    mean_maxClass = result_1$mean_maxClass 
    
    sd_minClass = result_1$sd_minClass
    sd_maxClass = result_1$sd_maxClass
    
    min_interval = result_1$min_interval 
    max_interval = result_1$max_interval 
    
    CHVP_test = result_1$CHVP_test 
    
    P_minClass = result_1$P_minClass 
    P_maxClass = result_1$P_maxClass
    
    P_minClass_mean = result_1$P_minClass_mean
    P_maxClass_mean = result_1$P_maxClass_mean
    
    P_minClass_permut = result_1$P_minClass_permut
    P_maxClass_permut = result_1$P_maxClass_permut
    
    
    Cup_chi2 = PCMR_correct(result_1,X_clump1$beta_hat_2,X_clump1$seb2)
    temp = function(c){
        return(sum((sort((pchisq(Cup_chi2,1*c,lower.tail = F))) - sort((seq(Cup_chi2)/length(Cup_chi2))))^2))
    }
    
    correct_factor = optimize(temp,lower=0,upper=10)$minimum
    CHVP_test_correct = pchisq(qchisq(result_1$CHVP_test,1,lower.tail = F),correct_factor,lower.tail = F) 
    
    print(c(CHVP_test,correct_factor, CHVP_test_correct))
    
    AIC_1 = result_1$AIC
    gamma_1 = str_c(result_1$gamma,collapse = "_")
    prop_1 = str_c(result_1$pi_gamma,collapse = "_")
    corTest_1 = result_1$Pvalue_corTest
  
    cat(paste(risk$Abbreviation,disease$Abbreviation,dim(X_clump)[1],P_1,eff_1,        
        mean_minClass,mean_maxClass,
        sd_minClass,sd_maxClass,
        min_interval,max_interval,
        CHVP_test, correct_factor, CHVP_test_correct,
        P_minClass,P_maxClass,
        P_minClass_mean,P_maxClass_mean,
        P_minClass_permut,P_maxClass_permut,
        AIC_1,gamma_1,prop_1,corTest_1,"\n",sep="\t"), file=saveFile,append=TRUE)
    print(paste0(S," end"))
}


Dir = "./Application/Psychiatric_disorders_analysis/trait/"
savePath = "./Application/Psychiatric_disorders_analysis/results/"
tests = read.csv(paste0(Dir,"common disease.csv")) 

saveFile = paste(savePath,"PCMR(n=2)_intact_boostrap_1000correct_log_05.txt",sep="/")
file.create(paste(saveFile,sep="/"))
cat(paste("risk","disease","IVs","P_1","eff_1",    
    "mean_minClass","mean_maxClass",
    "sd_minClass","sd_maxClass",
    "min_interval","max_interval",
    "CHVP_test","correct_factor", "CHVP_test_correct",
    "P_minClass","P_maxClass",
    "P_minClass_mean","P_maxClass_mean",
    "P_minClass_permut","P_minClass_permut",
    "AIC_1","gamma_1","prop_1","corTest_1\n",sep="\t"), file=saveFile)

            # cl <- makeCluster(10)
            # clusterExport(cl,"saveFile",envir = environment())
            # clusterExport(cl,"saveSvg",envir = environment())
            # clusterExport(cl,"boot",envir = environment())
            # clusterExport(cl,"PCMR_testCausal_bootstrap",envir = environment())
            # parLapply(cl,seq(dim(tests)[1]),sim_cause)
            # stopCluster(cl) # ?رռ?Ⱥ

for(i in seq(dim(tests)[1])){
    sim_cause(i)
}


