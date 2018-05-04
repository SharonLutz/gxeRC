# gxeRC
Power analyses for gene by environment interactions of both rare and common variants.

# Installation
```
install.packages("devtools") # devtools must be installed first

devtools::install_github("SharonLutz/gxeRC")
```

# Input
n is the number of subjects \
nSNP is the number of SNPS \
MAF is the minor allele frequency for the SNPS \
betaX is the genetic effect of each SNP \
betaI is the effect of interaction for each SNP \
zMu is the mean for the environmental effect \
zVar is the variance for the environmental effect \
yVar is the variance of the outcome Y \
nSim is the number of simulations \
alpha is the alpha level, default=0.05

# Example


```
library(gxeRC)
?gxeRC # For details on this function

gxeRC(n=5000)
```

# Output
For this example, we get the following matrix:

```
      lmX1  lmX2  lmX3 lmAll
[1,] 0.016 0.015 0.020 0.052
[2,] 0.095 0.030 0.017 0.154
[3,] 0.439 0.091 0.045 0.578
```
<img src="https://github.com/SharonLutz/gxeRC/blob/master/gxeRC.png" width="600">

# References
The power analysis used here is detailed in the following manuscript: <br/>
```
Lutz SM, Frederiksen B, Begum F, Cho MH, Hobbs B, McDonald ML, Parker
MM, DeMeo DL, Jiang L, Eringher M, Young K, Foreman MG, Kinney GL,
Make BJ, Lomas DA, Bakke P, Gulsvik A, Crapo JD, Silverman EK, Beaty
TH, Hokanson JE. (2018) Common and Rare Variants Analysis of Smoking
Related Traits Among Current and Former Smokers of European and
African Ancestry.  (Submitted).
```

