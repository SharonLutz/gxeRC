# gxeRC
The gxeRC R package computes power analyses for gene by environment interactions of both rare and common variants.

# Installation
```
install.packages("devtools") # devtools must be installed first

devtools::install_github("SharonLutz/gxeRC")
```

# Input
For n subjects, the number of SNPs X inputted by the user (input: nSNP) are genereated from binomial distributions with minor allele frequency specified by the user (input: MAF). The environmental factor Z is generated from a normal distribution with mean and variance inputted by the user. The outcome Y is generated from a normal distribution with mean as follows:

E\[Y\] = &beta;<sub>0</sub> + &beta;<sub>z</sub> Z + &sum;<sub>j</sub>  &beta;<sub>x</sub> X<sub>j</sub> + &sum;<sub>j</sub> &beta;<sub>I</sub> X<sub>j</sub>  Z  

for j=1,...,k where k= the number of SNPs.    

See the manpage for more detail regarding the input of the gxeRC function.

```
library(gxeRC)
?gxeRC # For details on this function
```

# Example
For 5,000 subjects, 3 SNPs X with MAF of 0.05, 0.01, and 0.005, respectively, and a normally distributed environmental factor Z, we genereated the power for the SNP by environment interaction on the outcome Y for each SNP independently and for the joint interaction using the following commands.

```
library(gxeRC)
gxeRC(n = 5000, nSNP = 3, MAF = c(0.05, 0.01, 0.005), betaX = c(0.25, 0.25, 0.25), 
betaI = c(0, 0.05, 0.1),zMu = 0, zVar = 1, yVar = 1, nSim = 1000, alpha = 0.05)
```

# Output
For this example, we get the following matrix and corresponding plot . We can see from the 

```
      lmX1  lmX2  lmX3 lmAll
[1,] 0.016 0.015 0.020 0.052
[2,] 0.095 0.030 0.017 0.154
[3,] 0.439 0.091 0.045 0.578
```
<img src="https://github.com/SharonLutz/gxeRC/blob/master/gxeRC.png" width="600">

# References
The power analysis used here was implemented in the following manuscript: <br/>

**Lutz SM**, Frederiksen B, Begum F, Cho MH, Hobbs B, McDonald ML, Parker
MM, DeMeo DL, Jiang L, Eringher M, Young K, Foreman MG, Kinney GL,
Make BJ, Lomas DA, Bakke P, Gulsvik A, Crapo JD, Silverman EK, Beaty
TH, Hokanson JE. (2018) Common and Rare Variants Analysis of Smoking
Related Traits Among Current and Former Smokers of European and
African Ancestry.  (Submitted).
