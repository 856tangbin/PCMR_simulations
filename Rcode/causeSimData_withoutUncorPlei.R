source("./Rcode/causeSimData.R")

Path = "./"
savePath = paste0(Path,"causeSimData_withoutUncorPlei/")

n1 = 40000
n2 = 40000
gamma = 0
q = 0
eta = 0

neffect2= 0 # setting the associated variant of exposure at zero in the absence of uncorrelated horizontal pleiotropy.

# In the absence of correlated horizontal pleiotropy 
for(tau in c(0,0.01,0.02,0.03,0.04,0.05)){
    for(omega in c(0)){
        gamma = sqrt(tau)
        eta = sqrt(omega)
        simData_ex(gamma,q,eta,n1=n1,n2=n2,neffect2=neffect2,seed=sum(as.integer(sapply(drop(stringr::str_split(paste0(gamma,q,eta,n1,n2),"",simplify = T)),charToRaw))))
    }
}

# In the presence of correlated horizontal pleiotropy
gamma = 0
for(q in c(0.1,0.2,0.3,0.4,0.5)){
    for(omega in c(0.01,0.02,0.03,0.04,0.05)){
        tau = gamma^2
        eta = sqrt(omega) 
        simData_ex(gamma,q,eta,n1=n1,n2=n2,neffect2=neffect2,seed=sum(as.integer(sapply(drop(stringr::str_split(paste0(gamma,q,eta,n1,n2),"",simplify = T)),charToRaw))))
    }
}
