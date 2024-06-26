## Generating GWAS summary statistics

We simulate GWAS summary statistics using the R package`causeSim`: https://github.com/jean997/causeSims. We added the function of setting random number seeds based on`causeSim`： `causeSims_0.1.0.tar.gz` to make the data repeatable. In this package, the LD structure dataset is estimated by 19,490 HapMap variants  on chromosome 19 in the CEU 1000 Genomes population, and this LD pattern is replicated 30 times, generating a genome-sized dataset of 584,700 variants. We apply the same parameters of the CAUSE method [[PMID: 32451458](https://pubmed.ncbi.nlm.nih.gov/32451458/)] in simulations, setting the heritability of exposure and outcome to 0.25, the number of susceptible variants of exposure and outcome to 1000, and the GWAS sample size of exposure and outcome to 40,000. Based on these parameters, we simulated the following scenarios:

| Shell script in this repository      | Scenarios                                                    |
| :----------------------------------- | :----------------------------------------------------------- |
| `1.causeSimData_withUncorPlei.sh`    | Simulate the null hypothesis ($`\gamma=0`$) in the presence of uncorrelated horizontal pleiotropy.<br />    $`q \in \{ 0\%,10\%,20\%,30\%,40\%,50\% \} `$ & $`\eta\in\{\sqrt{0.01},\sqrt{0.02},\sqrt{0.03},\sqrt{0.04},\sqrt{0.05}\}`$; |
| `2.causeSimData_withUncorPleiH1.sh`  | Simulate the alternative hypothesis in the presence of uncorrelated horizontal pleiotropy.<br />    1. In the absence of correlated horizontal pleiotropy. $`\gamma \in\{\sqrt{0.01},\sqrt{0.02},\sqrt{0.03}，\sqrt{0.04},\sqrt{0.05}\}`$;<br />    2. The correlated horizontal pleiotropic effect is in the same direction of the causal effect ($`\gamma=0.1`$). $`q\in \{10\%,20\%,30\%,40\%,50\%\}`$ & $`\omega^2+\tau^2=5\%`$; <br />    3. The correlated horizontal pleiotropic effect is in the opposite direction of the causal effect ($`\gamma=-0.1`$). $`q\in \{10\%,20\%,30\%,40\%,50\%\}`$ & $`\omega^2+\tau^2=5\%`$; |
| `3.causeSimData_withoutUncorPlei.sh` | Simulate the scenarios in the absence of uncorrelated horizontal pleiotropy ($`n_2=0`$).<br />    1. In the absence of correlated horizontal pleiotropy ($`\eta=0`$). $`\gamma \in\{\sqrt{0.01},\{\sqrt{0.00},\sqrt{0.02},\sqrt{0.03}，\sqrt{0.04},\sqrt{0.05}\}`$;<br />    2. In the presence of correlated horizontal pleiotropy under the null hypothesis ($`\gamma= 0`$). <br />    $`q\in \{0\%,10\%,20\%,30\%,40\%,50\%\}`$ & $`\eta\in\{\sqrt{0.01},\sqrt{0.02},\sqrt{0.03}，\sqrt{0.04},\sqrt{0.05}\}`$; |

The simulated GWAS summary statistics are saved in two folders, `./causeSimData` and `./causeSimData_withoutUncorPlei`. All the simulation data is so large that here we put only two examples in the folders. 

## Repeating the simulation analyses

We compared the performance of PCMR and other MR methods (IVW, Egger, Weighted median, Weighted mode, CAUSE and MRAID; MR-PRESSO) based on the generated GWAS summary statistics. The R packages of these methods compared in this study are shown in the following table with the version that we run.

|     Method      |                     R package                     |  Version   | Function  | parameters |
| :-------------: | :-----------------------------------------------: | :--------: | :-------: | :--------: |
|       IVW       |              MendelianRandomization               |   0.6.0    |  mr_ivw   |  default   |
|      Egger      |              MendelianRandomization               |   0.6.0    | mr_egger  |  default   |
| Weighted-median |              MendelianRandomization               |   0.6.0    | mr_median |  default   |
|   Median mode   |              MendelianRandomization               |   0.6.0    |  mr_mbe   |  default   |
|      CAUSE      |     [cause](https://github.com/jean997/cause)     | 1.2.0.0335 |   cause   |  default   |
|      MRAID      | [MRAID](https://github.com/yuanzhongshang/MRAID)  |    1.0     |   MRAID   |  default   |
|    MR-PRESSO    | [MRPRESSO](https://github.com/rondolab/MR-PRESSO) |    1.0     | mr_presso |  default   |
|      PCMR       |    [PCMR](https://github.com/856tangbin/PCMR)     |   0.1.0    |   PCMR    |  default   |

Paths to the executables of the methods are in `./MRmethod`, and all results are provided in the files,`./results` and `./results_withoutUncorPlei`. 

## Real data analyses

The codes for analyses of each method for common diseases and psychiatric disorders, and the results of the analyses are placed in the file`./Application`. For PCMR, we also plotted the clustering results under various trait pairs, at `./Application/Common_disease_analysis/results/PCMR(n=2)_intact` and  `./Application/Psychiatric_disorders_analysis/results/PCMR(n=2)_intact`. 

