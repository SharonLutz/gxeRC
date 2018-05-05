gxeRC <-
function(n=5000,nSNP=3,MAF=c(0.05,0.01,0.005),betaX=c(0.25,0.25,0.25),betaI=c(0,0.05,0.1),
                zMu=0,zVar=1,yVar=1,nSim=1000,alpha=0.05,plot.name="gxeRC.pdf"){ 
  
  ####################################
  # input parameters
  ####################################
  # n is the number of subjects
  # nSNP= the number of SNPS
  # MAF= minor allele frequency for the SNPS
  # betaX = genetic effect of each SNP
  # betaI= effect of interaction for each SNP
  # zMu is the mean for the environmental effect
  # zVar is the variance for the environmental effect
  
  # nSim is the number of simulations
  # alpha is the alpha level, default=0.05
  
  ####################################
  # Error Checks
  ####################################
  
  # Check nSNP = length(MAF)==length(betaX)==length(betaI)
  if(nSNP!=length(MAF)){stop("Error: nSNP must equal length(MAF).")}
  if(nSNP!=length(betaX)){stop("Error: nSNP must equal length(betaX).")}
  if(nSNP!=length(betaI)){stop("Error: nSNP must equal length(betaI).")}
  
  # Check n, nSNP and nSim are integers
  if(floor(n)!=ceiling(n)){stop("Error: n must be an integer.")}
  if(floor(nSNP)!=ceiling(nSNP)){stop("Error: nSNP must be an integer.")}
  if(floor(nSim)!=ceiling(nSim)){stop("Error: nSim must be an integer.")}
  
  # Check n, nSNP and nSim are greater than 0
  if(!(n>0)){stop("Error: n must be greater than 0.")}
  if(!(nSNP>0)){stop("Error: nSNP must be greater than 0.")}
  if(!(nSim>0)){stop("Error: nSim must be greater than 0.")}
  
  # Check zVar > 0 and yVar > 0
  if(!(zVar>0)){stop("Error: zVar must be greater than 0.")}
  if(!(yVar>0)){stop("Error: yVar must be greater than 0.")}
  
  # Check length(zVar)==1 length(zMu)==1 & length(yVar==1)
  if(length(zVar)!=1){stop("Error: zVar must be of length 1")}
  if(length(zMu)!=1){stop("Error: zMu must be of length 1")}
  if(length(yVar)!=1){stop("Error: yVar must be of length 1")}
  
  # Check alpha>0 & alpha<1
  if(alpha<0 | alpha>1){stop("Error: alpha must be between 0 and 1.")}
  
  ####################################
  # Store Results
  #################################### 
  
  rejectH0<-matrix(0,nrow=length(betaI),ncol=(nSNP+1))
  colnames(rejectH0)<-c(paste("lmX",1:nSNP,sep=""),"lmAll")
  
  ####################################
  # Run Simulations
  #################################### 
  
  for(GLOBALVAR in 1:nSim){
    set.seed(GLOBALVAR)
    if(floor(GLOBALVAR/100)==ceiling(GLOBALVAR/100)){print(GLOBALVAR)}
    
    ####################################
    # Generate Data
    ####################################
    
    # CYCLE through values of betaI
    for(bb in 1:length(betaI)){
      betaIv<-betaI[bb]
      
      ####################################
      # simulate data
      ####################################
      
      # generate the matrix of SNPs
      X<-matrix(0,nrow=n,ncol=nSNP)
      
      errorFound <- F
      for(xx in 1:nSNP){
        X[,xx]<-rbinom(n,2,MAF[xx]) 
        # Check X is not all zero
        if(mean(X[,xx])==0|mean(X[,xx])==2){
          problemSNP<- xx
          errorFound <- T
          break
          } # let user know they need to increase n or MAF because there is no variability in SNP xx <- give what xx is don't give xx as an index
      }
      if(errorFound){
        errormessage <- paste("Error: Increase n or MAF because there is no variability in SNP ",problemSNP,sep = "")
        stop(errormessage)}
      
      # generate the environment Z
      z<- rnorm(n,zMu,zVar) 
      zz<-matrix(0,nrow=n,ncol=1)
      zz[,1]<-z
      
      # generate the outcome Y
      mainEffects<- X%*%betaX
      intEffects<- (X%*%rep(betaIv,nSNP))*zz
      yMu<- mainEffects+ intEffects
      y<-rnorm(n,yMu,yVar)
      
      ####################################
      # linear regression for interaction
      ####################################
      modelA<-lm(y~z+X+X*z)
      modelAA<-summary(modelA)$coef
      
      if(nrow(modelAA)<(nSNP+2+nSNP)){stop("Error: Increase n or MAF because there is not enough variability")}
      
      nRow<-nSNP+2
      
      for(rr in 1:nSNP){
        if(modelAA[(nRow+rr),4]<(alpha/nSNP)){rejectH0[bb,rr]<-rejectH0[bb,rr]+1}
      }
      
      modelR<-lm(y~z+X)
      if(anova(modelA,modelR)$P[2]<alpha){rejectH0[bb,"lmAll"]<-rejectH0[bb,"lmAll"]+1}
      
      ####################################
      # Compile results 
      ####################################
    }#end betaI loop
  }#end globalvar
  
  rejectMat<-rejectH0/nSim
  
  ####################################
  # Create plot
  #################################### 
  
  nn<-ncol(rejectMat)
  
  pdf(plot.name)
  plot(-1,-1,xlim=c(min(betaI),max(betaI)),ylim=c(0,1),xlab="betaI",ylab="",main="")
  for(pp in 1:nn){
    lines(betaI,rejectMat[,pp],pch=pp,col=pp,type="b")
  }
  legend("topleft",c(paste("SNP",1:nSNP,": MAF=",c(MAF),sep=""),"All SNPs"),col=c(1:nn),pch=(1:nn),lwd=1)
  dev.off()
  
  ####################################
  # End function
  ####################################
  
  list(rejectMat)}
