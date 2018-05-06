# gxeRC
Power analyses for gene by environment interactions of both rare and common variants.

# Installation
```
install.packages("devtools") # devtools must be installed first

devtools::install_github("SharonLutz/gxeRC")
```

# Input

nSNP is the number of SNPs generated from a binomial distribution for n subjects (input n) for a given minor allele frequency (input vector MAF).

# Example


```
library(gxeRC)
?gxeRC # For details on this function

gxeRC(n=5000)
```

# Output
For this example, we get the following matrix and corresponding plot:

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


