\name{gxeRC}
\alias{gxeRC}
\title{gxeRC
}
\description{
A function to  examine gene by environment interactions of both rare and common variants.
}
\usage{
gxeRC(n = 5000, nSNP = 3, MAF = c(0.05, 0.01, 0.005), betaX = c(0.25, 0.25, 0.25), betaI = c(0, 0.05, 0.1), 
zMu = 0, zVar = 1, yVar = 1, nSim = 1000, alpha = 0.05, plot.name = "gxeRC.pdf")
}

\arguments{
  \item{n}{
is the number of subjects
}
  \item{nSNP}{
is the number of SNPS
}
  \item{MAF}{
is the minor allele frequency for the SNPS
}
  \item{betaX}{
is the genetic effect of each SNP
}
  \item{betaI}{
is the effect of interaction for each SNP
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
a matrix of the proportion of simulations the null hypothesis was rejected for each SNP, and all SNPs together
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
gxeRC(n=5000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gene }% use one of  RShowDoc("KEYWORDS")
\keyword{ interaction }% __ONLY ONE__ keyword per line
