source("./Rcode/causeSimData.R")

n1 = 40000
n2 = 40000


# H1 in absence of correlated pleiotropy
for(tau in c(0.01,0.02,0.03,0.04,0.05)){
    for(q in c(0)){
        for(omega in c(0)){
            gamma = sqrt(tau)
            eta = 0
            simData_ex(gamma,q,eta,n1,n2,seed=sum(as.integer(sapply(drop(stringr::str_split(paste0(gamma,q,eta,n1,n2),"",simplify = T)),charToRaw))))
        }
    }
}

# H1 in presence of correlated pleiotropy. The correlated horizontal pleiotropic effect is in the same direction of the causal effect


gamma = -0.1
q = 0
eta = 0
simData_ex(gamma,q,eta,n1=n1,n2=n2,seed=sum(as.integer(sapply(drop(stringr::str_split(paste0(gamma,q,eta,n1,n2),"",simplify = T)),charToRaw))))
for(q in c(0.1,0.2,0.3,0.4,0.5)){
    omega = 0.05 - gamma**2
    eta = sqrt(omega/q)
    simData_ex(gamma,q,eta,n1=n1,n2=n2,seed=sum(as.integer(sapply(drop(stringr::str_split(paste0(gamma,q,eta,n1,n2),"",simplify = T)),charToRaw))))
}

# H1 in presence of correlated pleiotropy. The correlated horizontal pleiotropic effect is in the opposite direction of the causal effect
gamma = 0.1
# simData_ex(gamma,q,eta,n1=n1,n2=n2,seed=sum(as.integer(sapply(drop(stringr::str_split(paste0(gamma,q,eta,n1,n2),"",simplify = T)),charToRaw)))) # already running
for(q in c(0.1,0.2,0.3,0.4,0.5)){
    omega = 0.05 - gamma**2
    eta = sqrt(omega/q)
    simData_ex(gamma,q,eta,n1=n1,n2=n2,seed=sum(as.integer(sapply(drop(stringr::str_split(paste0(gamma,q,eta,n1,n2),"",simplify = T)),charToRaw))))
}