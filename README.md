# PCMR_simulations

## Generating GWAS summary statistics

我们使用causeSim模拟GWAS summary statistics，用于比较不同MR方法的效果：

We simulate summary statistics using R package `causeSim`: https://github.com/jean997/causeSims. In order to make the data repeatable, we added the function of setting random number seeds on the basis of `causeSim`： `causeSims_0.1.0.tar.gz`. In this package, the LD structure dataset is estimated by 19,490 HapMap variants  on chromosome 19 in the CEU 1000 Genomes population, and this LD pattern is replicated 30 times, generating a genome-sized dataset of 584,700 variants. In simulations, we apply the same parameters of the CAUSE method [[PMID: 32451458](https://pubmed.ncbi.nlm.nih.gov/32451458/)], setting the heritability of exposure and outcome to 0.25, the number of susceptible variants of exposure and outcome to 1000, and the GWAS sample size of exposure and outcome to 40,000. 

| Shell script in this repository      | scenarios                                                    |
| :----------------------------------- | :----------------------------------------------------------- |
| `1.causeSimData_withUncorPlei.sh`    | Simulate the null hypothesis ($\gamma=0$) in the presence of uncorrelated horizontal pleiotropy.<br />$q \in \{ 0\%,10\%,20\%,30\%,40\%,50\%\}$ <br />& $\eta\in\{\sqrt{0.01},\sqrt{0.02},\sqrt{0.03},\sqrt{0.04},\sqrt{0.05}\}$; |
| `2.causeSimData_withUncorPleiH1.sh`  | Simulate the alternative hypothesis in the presence of uncorrelated horizontal pleiotropy.<br />1. In the absence of correlated horizontal pleiotropy. $\gamma \in\{\sqrt{0.01},\sqrt{0.02},\sqrt{0.03}，\sqrt{0.04},\sqrt{0.05}\}$;<br />2. The correlated horizontal pleiotropic effect is in the same direction of the causal effect ($\gamma=0.1$). $q\in \{10\%,20\%,30\%,40\%,50\%\}$ & $\omega^2+\tau^2=5\%$; <br />3. The correlated horizontal pleiotropic effect is in the opposite direction of the causal effect ($\gamma=-0.1$). $q\in \{10\%,20\%,30\%,40\%,50\%\}$ & $\omega^2+\tau^2=5\%$; |
| `3.causeSimData_withoutUncorPlei.sh` | Simulate the scenarios in the absence of uncorrelated horizontal pleiotropy ($n_2=0$).<br />1. In the absence of correlated horizontal pleiotropy ($\eta=0$). $\gamma \in\{\sqrt{0.01},\{\sqrt{0.00},\sqrt{0.02},\sqrt{0.03}，\sqrt{0.04},\sqrt{0.05}\}$;<br />2. In the presence of correlated horizontal pleiotropy under the null hypothesis ($\gamma= 0$). <br />$q\in \{0\%,10\%,20\%,30\%,40\%,50\%\}$ & $\eta\in\{\sqrt{0.01},\sqrt{0.02},\sqrt{0.03}，\sqrt{0.04},\sqrt{0.05}\}$; |

The simulated GWAS summary statistics are saved in two folders, `./causeSimData` and `./causeSimData_withoutUncorPlei`. All the simulation data is so large that here we put only two examples in the folders. 

$$
q \in \{0\%,10\%,20\%,30\%,40\%,50\%\}
$$


## 2. 

