library(cause)
library(parallel) # ???߳?


sim_cause = function(i){
    library(cause)
    library(readr)
    library(dplyr)
    
    Dir = "./Psychiatric_disorders_analysis/trait/"
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
    
    # CAUSE process--------------------------
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
                    

    # Calculate nuisance parameters
    set.seed(100)
    varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
    params <- est_cause_params(X, varlist)
    
    r2_thresh = 0.1
    pval_thresh = 1e-3
    
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
    top_vars <- X_clump$rsid
    res <- cause(X=X, variants = top_vars, param_ests = params)
    result = summary(res,ci_size = 0.95)

    gamma_median_causal = result$quants[[2]][1,1]
    q_median_causal = result$quants[[2]][1,3]
    eta_median_causal = result$quants[[2]][1,2]
    
    q_median_sharing = result$quants[[1]][1,3]
    eta_median_sharing = result$quants[[1]][1,2]
    p_cause = result$p
    
    rho = params$rho
    
    cat(paste(risk$Abbreviation,disease$Abbreviation,dim(X_clump)[1],gamma_median_causal,	q_median_causal,	eta_median_causal,	q_median_sharing,	eta_median_sharing,	p_cause,	rho,"\n",sep="\t"), file=saveFile,append=TRUE)
}


Dir = "./Psychiatric_disorders_analysis/trait/"
savePath = "./Psychiatric_disorders_analysis/results/"
tests = read.csv(paste0(Dir,"common disease.csv")) 

saveFile = paste(savePath,"causeThred.txt",sep="/")
file.create(paste(saveFile,sep="/"))
cat(paste("risk","disease","IVs","gamma_median_causal","q_median_causal","eta_median_causal","q_median_sharing","eta_median_sharing","p_cause","rho\n",sep="\t"), file=saveFile)

cl <- makeCluster(10)
clusterExport(cl,"saveFile",envir = environment())
parLapply(cl,seq(dim(tests)[1]),sim_cause)

stopCluster(cl) # ?رռ?Ⱥ
