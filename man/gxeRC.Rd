\name{gxeRC}
\alias{gxeRC}
\title{gxeRC
}
\description{
A function to run a power analysis for gene by environment interactions of both rare and common variants.
}
\usage{
gxeRC(n = 5000, nSNP = 3, MAF = c(0.05, 0.01, 0.005), betaX = c(0.25, 0.25, 0.25), betaI = c(0, 0.05, 0.1), zMu = 0, zVar = 1, yVar = 1, nSim = 1000, alpha = 0.05, plot.name = "gxeRC.pdf")
}

\arguments{
  \item{n}{
is the number of subjects
}
  \item{nSNP}{
the number of SNPS
}
  \item{MAF}{
minor allele frequency for the SNPS
}
  \item{betaX}{
genetic effect of each SNP
}
  \item{betaI}{
effect of interaction for each SNP
}
  \item{zMu}{
is the mean for the environmental effect
}
  \item{zVar}{
is the variance for the environmental effect
}
  \item{yVar}{
is the variance for the outcome Y
}
  \item{nSim}{
is the number of simulations
}
  \item{alpha}{
is the alpha level, default=0.05
}
  \item{plot.name}{
is the name of the corresponding plot
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a matrix of the power of each SNP, and all SNPs together
}
\references{
Lutz SM, Frederiksen B, Begum F, Cho MH, Hobbs B, McDonald ML, Parker
MM, DeMeo DL, Jiang L, Eringher M, Young K, Foreman MG, Kinney GL,
Make BJ, Lomas DA, Bakke P, Gulsvik A, Crapo JD, Silverman EK, Beaty
TH, Hokanson JE. (2018) Common and Rare Variants Analysis of Smoking
Related Traits Among Current and Former Smokers of European and
African Ancestry.  (Submitted).
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n = 5000, nSNP = 3, MAF = c(0.05, 0.01, 0.005), betaX = c(0.25, 
    0.25, 0.25), betaI = c(0, 0.05, 0.1), zMu = 0, zVar = 1, 
    yVar = 1, nSim = 1000, alpha = 0.05, plot.name = "gxeRC.pdf") 
{
    if (nSNP != length(MAF)) {
        stop("Error: nSNP must equal length(MAF).")
    }
    if (nSNP != length(betaX)) {
        stop("Error: nSNP must equal length(betaX).")
    }
    if (nSNP != length(betaI)) {
        stop("Error: nSNP must equal length(betaI).")
    }
    if (length(n) > 1 | !(n > 0)) {
        stop("Error: n must be an integer greater than 0.")
    }
    if (length(nSNP) > 1 | !(nSNP > 0)) {
        stop("Error: nSNP must be an integer greater than 0.")
    }
    if (length(nSim) > 1 | !(nSim > 0)) {
        stop("Error: nSim must be an integer greater than 0.")
    }
    if (!(zVar > 0)) {
        stop("Error: zVar must be greater than 0.")
    }
    if (!(yVar > 0)) {
        stop("Error: yVar must be greater than 0.")
    }
    if (length(zVar) != 1) {
        stop("Error: zVar must be of length 1")
    }
    if (length(zMu) != 1) {
        stop("Error: zMu must be of length 1")
    }
    if (length(yVar) != 1) {
        stop("Error: yVar must be of length 1")
    }
    if (alpha < 0 | alpha > 1) {
        stop("Error: alpha must be between 0 and 1.")
    }
    rejectH0 <- matrix(0, nrow = length(betaI), ncol = (nSNP + 
        1))
    colnames(rejectH0) <- c(paste("lmX", 1:nSNP, sep = ""), "lmAll")
    for (GLOBALVAR in 1:nSim) {
        set.seed(GLOBALVAR)
        if (floor(GLOBALVAR/100) == ceiling(GLOBALVAR/100)) {
            print(GLOBALVAR)
        }
        for (bb in 1:length(betaI)) {
            betaIv <- betaI[bb]
            X <- matrix(0, nrow = n, ncol = nSNP)
            errorFound <- F
            for (xx in 1:nSNP) {
                X[, xx] <- rbinom(n, 2, MAF[xx])
                if (mean(X[, xx]) == 0 | mean(X[, xx]) == 2) {
                  problemSNP <- xx
                  errorFound <- T
                  break
                }
            }
            if (errorFound) {
                errormessage <- paste("Error: Increase n or MAF because there is no variability in SNP ", 
                  problemSNP, sep = "")
                stop(errormessage)
            }
            z <- rnorm(n, zMu, zVar)
            zz <- matrix(0, nrow = n, ncol = 1)
            zz[, 1] <- z
            mainEffects <- X \%*\% betaX
            intEffects <- (X \%*\% rep(betaIv, nSNP)) * zz
            yMu <- mainEffects + intEffects
            y <- rnorm(n, yMu, yVar)
            modelA <- lm(y ~ z + X + X * z)
            modelAA <- summary(modelA)$coef
            if (nrow(modelAA) < (nSNP + 2 + nSNP)) {
                stop("Error: Increase n or MAF because there is not enough variability")
            }
            nRow <- nSNP + 2
            for (rr in 1:nSNP) {
                if (modelAA[(nRow + rr), 4] < (alpha/nSNP)) {
                  rejectH0[bb, rr] <- rejectH0[bb, rr] + 1
                }
            }
            modelR <- lm(y ~ z + X)
            if (anova(modelA, modelR)$P[2] < alpha) {
                rejectH0[bb, "lmAll"] <- rejectH0[bb, "lmAll"] + 
                  1
            }
        }
    }
    rejectMat <- rejectH0/nSim
    nn <- ncol(rejectMat)
    pdf(plot.name)
    plot(-1, -1, xlim = c(min(betaI), max(betaI)), ylim = c(0, 
        1), xlab = "betaI", ylab = "", main = "")
    for (pp in 1:nn) {
        lines(betaI, rejectMat[, pp], pch = pp, col = pp, type = "b")
    }
    legend("topleft", c(paste("SNP", 1:nSNP, ": MAF=", c(MAF), 
        sep = ""), "All SNPs"), col = c(1:nn), pch = (1:nn), 
        lwd = 1)
    dev.off()
    list(rejectMat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gene }% use one of  RShowDoc("KEYWORDS")
\keyword{ interaction }% __ONLY ONE__ keyword per line
