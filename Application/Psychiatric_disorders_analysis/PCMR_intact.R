library(cause)
library(parallel) # ???߳?

saveSvg = function(result,S,model){
    svg(
        filename = paste0("./Application/Psychiatric_disorders_analysis/results/PCMR(n=2)_intact/",S[1],"_",S[2],"model",model,".svg"), # ?ļ?????
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
    library(weights)
    
    library(stringr)
    
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
    
    #if(risk$Abbreviation == "height"){
    #   return()
    #}
    
    
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
    num_gamma = 2
    
    # model 1 : gSigma2 T ; sigma2 T
    result_1 = PCMR(X_clump$beta_hat_1, X_clump$seb1,
                    X_clump$beta_hat_2,X_clump$seb2,num_gamma = num_gamma,model="1",
                    sigma2 = init$sigma2,rho=init$rho,isIntact = T)  
    result_1 = PCMR_testCorPlei(result_1)
    saveSvg(result_1,S,1)

    # model 2 : gSigma2 F ; sigma2 T
    result_2 =  PCMR(X_clump$beta_hat_1, X_clump$seb1,
                    X_clump$beta_hat_2,X_clump$seb2,num_gamma = num_gamma,model="2",
                    sigma2 = init$sigma2,rho=init$rho,isIntact = T)
    result_2 = PCMR_testCorPlei(result_2)
    saveSvg(result_2,S,2)
    
    # model 3 : gSigma2 T ; sigma2 F
    result_3 =  PCMR(X_clump$beta_hat_1, X_clump$seb1,
                    X_clump$beta_hat_2,X_clump$seb2,num_gamma = num_gamma,model="3",
                    sigma2 = init$sigma2,rho=init$rho,isIntact = T)  
    result_3 = PCMR_testCorPlei(result_3)
    saveSvg(result_3,S,3)
    
    # model 4 : gSigma2 F ; sigma2 F
    result_4 =  PCMR(X_clump$beta_hat_1, X_clump$seb1,
                    X_clump$beta_hat_2,X_clump$seb2,num_gamma = num_gamma,model="4",
                    sigma2 = init$sigma2,rho=init$rho,isIntact = T)
    result_4 = PCMR_testCorPlei(result_4)
    saveSvg(result_4,S,4)
    
    # save results
    P_1 = result_1$Pvalue
    eff_1 = result_1$effect
    AIC_1 = result_1$AIC
    gamma_1 = str_c(result_1$gamma,collapse = "_")
    prop_1 = str_c(result_1$pi_gamma,collapse = "_")
    corTest_1 = result_1$Pvalue_corTest
    
    b_1 = str_c(result_1$b,collapse = "_")
    sigma2_1 = str_c(result_1$sigma2,collapse = "_")
    b_prop_1 = str_c(result_1$pi_b,collapse = "_")
    
    gSigma2_1 = str_c(result_1$gSigma2,collapse = "_")
    
  
    P_2 = result_2$Pvalue
    eff_2 = result_2$effect
    AIC_2 = result_2$AIC
    gamma_2 = str_c(result_2$gamma,collapse = "_")
    prop_2 = str_c(result_2$pi_gamma,collapse = "_")
    corTest_2 = result_2$Pvalue_corTest
        
    P_3 = result_3$Pvalue
    eff_3 = result_3$effect
    AIC_3 = result_3$AIC
    gamma_3 = str_c(result_3$gamma,collapse = "_")
    prop_3 = str_c(result_3$pi_gamma,collapse = "_")
    corTest_3 = result_3$Pvalue_corTest
    
    P_4 = result_4$Pvalue
    eff_4 = result_4$effect
    AIC_4 = result_4$AIC
    gamma_4 = str_c(result_4$gamma,collapse = "_")
    prop_4 = str_c(result_4$pi_gamma,collapse = "_")
    corTest_4 = result_4$Pvalue_corTest
    
    # cat(paste(risk$Abbreviation,disease$Abbreviation,dim(X_clump)[1],P_1,eff_1,AIC_1,gamma_1,prop_1,corTest_1,P_2,eff_2,AIC_2,gamma_2,prop_2,corTest_2,P_3,eff_3,AIC_3,gamma_3,prop_3,corTest_3,P_4,eff_4,AIC_4,gamma_4,prop_4,corTest_4,"\n",sep="\t"), file=saveFile,append=TRUE)
    cat(paste(risk$Abbreviation,disease$Abbreviation,dim(X_clump)[1],P_1,eff_1,AIC_1,gamma_1,gSigma2_1,prop_1,b_1,sigma2_1,b_prop_1,corTest_1,P_2,eff_2,AIC_2,gamma_2,prop_2,corTest_2,P_3,eff_3,AIC_3,gamma_3,prop_3,corTest_3,P_4,eff_4,AIC_4,gamma_4,prop_4,corTest_4,"\n",sep="\t"), file=saveFile,append=TRUE)
    print(paste0(S," end"))
    
    
    # ???????? ##############################
    result = result_1
    class1 = (rowSums(result$W_ind,dims = 2)[,1] == 1)
    class2 = (rowSums(result$W_ind,dims = 2)[,2] == 1)
    
    # class1 has the argmax q
    if(result$pi_gamma[1] < result$pi_gamma[2]){
        t = class1
        class1 = class2
        class2 = t
    }
    
    
    # # Calculating LD to determine which class each snp belong to. #################
    # 
    # variants = X[X$p1 < pval_thresh,"snp"]
    # 
    # # ֱ??ʹ?? pruning variants  
    # # significant_class1 = X_clump$rsid[class1]
    # # significant_class2 = X_clump$rsid[class2]
    # 
    # # ʹ????????ѡ????ͬ???Ĺ??߱?��
    # # ?ж? variants ?ĳ??ȣ?????̫??????Ҫ?ֶμ??㣬?????Ϊ5000һ??
    # Len = length(variants)
    # 
    # sep = 1000
    # blocks = Len %/% sep + 1
    # significant_class1 = c()
    # significant_class2 = c()
    # 
    # for(j in seq(blocks)){
    #     start = (j-1) * sep + 1
    #     
    #     if(j == blocks){
    #         end = Len
    #     }else{
    #         end = j * sep
    #     }
    #     
    #     print(c(start,end))
    #     
    #     
    #     variants_block = c(variants[start:end],X_clump$rsid)
    #     matrix_LD = ieugwasr::ld_matrix(variants_block ,
    #                                     with_alleles = F,
    #                                     plink_bin = "/home/tb/tools/GWAS//plink",
    #                                     bfile="/home/tb/ref/1kg.v3//EUR")
    #     variants_block = colnames(matrix_LD) # ??һЩ?޷?????LD?ĸ??޳???
    #     
    #     LD_class1 = matrix_LD[,X_clump$rsid[class1]]
    #     LD_class2 = matrix_LD[,X_clump$rsid[class2]]
    #     
    #     # ???????Ե?????ϵ????С
    #     distance_class1 = apply(LD_class1 ** 2,1,max)
    #     distance_class2 = apply(LD_class2 ** 2,1,max)
    #     
    #     
    #     significant_class1 = c(significant_class1,variants_block[distance_class1 >= distance_class2])
    #     significant_class2 = c(significant_class2,variants_block[distance_class1 <= distance_class2])
    # }
    # 
    # c1 = risk_data[risk_data[[risk$rsid]] %in% significant_class1,c(risk$rsid,risk$CHR,risk$BP,risk$P.value)]
    # colnames(c1) = c("rsid","chromosome","base_pair_location","P-value")
    # c1$class = 1
    # 
    # c2 = risk_data[risk_data[[risk$rsid]] %in% significant_class2,c(risk$rsid,risk$CHR,risk$BP,risk$P.value)]
    # colnames(c2) = c("rsid","chromosome","base_pair_location","P-value")
    # c2$class = 2
    # 
    # if(any(c1$chromosome == 0) | any(c2$chromosome == 0)){
    #     c1$chromosome = c1$chromosome + 1
    #     c2$chromosome = c2$chromosome + 1
    # }
    # 
    # write.table(c1[,seq(2,4)],paste0("./Application/enrichAnalysis_class/traits//",i,"_",risk$Abbreviation,"_",disease$Abbreviation,"_class1.txt"),row.names = FALSE,sep=",",quote = FALSE)
    # write.table(c2[,seq(2,4)],paste0("./Application/enrichAnalysis_class/traits//",i,"_",risk$Abbreviation,"_",disease$Abbreviation,"_class2.txt"),row.names = FALSE,sep=",",quote = FALSE)
    # 
    # c3 = rbind(c1,c2)
    # c3[c3$rsid %in% X_clump$rsid,"isPruning"] = 1
    # c3$risk = risk$Abbreviation
    # c3$disease = disease$Abbreviation
    # 
    # write.csv(c3,paste0("./Application/enrichAnalysis_class/traits//",i,"_",risk$Abbreviation,"_",disease$Abbreviation,".csv"),row.names = FALSE,quote = FALSE)
}


Dir = "./Application/Psychiatric_disorders_analysis/trait/"
savePath = ".lication/Psychiatric_disorders_analysis/results/"
tests = read.csv(paste0(Dir,"common disease.csv")) 

saveFile = paste(savePath,"PCMR(n=2)_intact.txt",sep="/")
file.create(paste(saveFile,sep="/"))
# cat(paste("risk","disease","IVs","P_1","eff_1","AIC_1","gamma_1","prop_1","corTest_1","P_2","eff_2","AIC_2","gamma_2","prop_2","corTest_2","P_3","eff_3","AIC_3","gamma_3","prop_3","corTest_3","P_4","eff_4","AIC_4","gamma_4","prop_4","corTest_4\n",sep="\t"), file=saveFile)

cat(paste("risk","disease","IVs","P_1","eff_1","AIC_1","gamma_1","gSigma2_1","prop_1","b_1","sigma2_1","b_prop_1","corTest_1","P_2","eff_2","AIC_2","gamma_2","prop_2","corTest_2","P_3","eff_3","AIC_3","gamma_3","prop_3","corTest_3","P_4","eff_4","AIC_4","gamma_4","prop_4","corTest_4\n",sep="\t"), file=saveFile)

cl <- makeCluster(48)
clusterExport(cl,"saveFile",envir = environment())
clusterExport(cl,"saveSvg",envir = environment())
parLapply(cl,seq(dim(tests)[1]),sim_cause)

stopCluster(cl) # ?رռ?Ⱥ
