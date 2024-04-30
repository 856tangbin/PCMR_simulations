library(cause)
library(parallel) # ???߳?


sim_cause = function(i){
    library(cause)
    library(readr)
    library(dplyr)
    library(MendelianRandomization)
    
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
    
    # model fitting
    # MendelianRandomization -------------------------
    mr_data = mr_input(bx = X_clump$beta_hat_1, bxse = X_clump$seb1,
                       by = X_clump$beta_hat_2, byse = X_clump$seb2)
    
    result = mr_mbe(mr_data)
    P_mbe = result$Pvalue
    eff_mbe = result$Estimate
    se_mbe = result$StdError
    
    result = mr_egger(mr_data)
    P_egger = result$Pvalue.Est
    eff_egger = result$Estimate
    se_egger = result$StdError.Est
    
    result = mr_ivw(mr_data)
    P_ivw = result$Pvalue
    eff_ivw = result$Estimate
    se_ivw = result$StdError
    
    result = mr_median(mr_data)
    P_weightedMedian = result$Pvalue
    eff_weightedMedian = result$Estimate
    se_weightedMedian = result$StdError
    
    # save results
    cat(paste(risk$Abbreviation,disease$Abbreviation,dim(X_clump)[1],P_mbe,eff_mbe,se_mbe,P_egger,eff_egger,se_egger,P_ivw,eff_ivw,se_ivw,P_weightedMedian,eff_weightedMedian,se_weightedMedian,"\n",sep="\t"), file=saveFile,append=TRUE)
}


Dir = "./Application/Psychiatric_disorders_analysis/trait/"
savePath = "./Application/Psychiatric_disorders_analysis/results/"
tests = read.csv(paste0(Dir,"common disease.csv")) 

saveFile = paste(savePath,"mendelianRandomization.txt",sep="/")
file.create(paste(saveFile,sep="/"))
cat(paste("risk","disease","IVs","P_mbe","eff_mbe","se_mbe","P_egger","eff_egger","se_egger","P_ivw","eff_ivw","se_ivw","P_weightedMedian","eff_weightedMedian","se_weightedMedian\n",sep="\t"), file=saveFile)

cl <- makeCluster(48)
clusterExport(cl,"saveFile",envir = environment())
parLapply(cl,seq(dim(tests)[1]),sim_cause)

stopCluster(cl) # ?رռ?Ⱥ
