source("./Rcode/causeSimData.R")

n1 = 40000 # GWAS sample size of exposure
n2 = 40000 # GWAS sample size of outcome
gamma = 0 # causal effect
q = 0 # the proportion of correlated horizontal pleiotropic effect
eta = 0 # correlated horizontal pleiotropic effect 

# H0 without correlated pleiotropy ###################################################
simData_ex(gamma,q,eta,n1=n1,n2=n2,seed=sum(as.integer(sapply(drop(stringr::str_split(paste0(gamma,q,eta,n1,n2),"",simplify = T)),charToRaw))))

# H0 with correlated pleiotropy ##############################################################
for(gamma in c(0)){
    for(q in c(0.1,0.2,0.3,0.4,0.5)){
        for(omega in c(0.01,0.02,0.03,0.04,0.05)){
            eta = sqrt(omega)
            simData_ex(gamma,q,eta,n1=n1,n2=n2,seed=sum(as.integer(sapply(drop(stringr::str_split(paste0(gamma,q,eta,n1,n2),"",simplify = T)),charToRaw))))
        }
    }
}

