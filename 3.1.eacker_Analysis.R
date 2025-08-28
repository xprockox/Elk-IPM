### Elk Integrated Population Model Adapted from Eacker et al. 2017
### Last updated: July 10, 2025
### Contact: xprockox@gmail.com

############################################################################################
### Load packages
library(nimble)
library(MCMCvis)
library(dplyr)

#####################################################################################################
# Data, intial values, parameters to save, and other code to execute IPM in JAGS using
# R2jags package in Program R for East Fork (EF) elk population in Bitterroot Valley
#
#
#####################################################################################################
#####################################################################################################
# Data for EF base IPM

#####################################################################################################
# summer calf survival data (both sexes)

# here, each row represents a year.
# since this is calf survival, there is only one year of monitoring per individual. therefore, consecutive
# years in the same column do not correspond to the same individual. instead, each column  
# represents a "slot" that an individual could occupy. for the first year,
# only 44 individuals were monitored. this is why there are only 44 complete columns. 
# the next year (year 2) included 53 monitored individuals. this is why there are only 53 complete columns. 



#exit times (when individuals die (days))
S.Y.S<-structure(.Data=c(NA, 82, NA, 7, NA, NA, NA, NA, NA, NA, 10, 47, 69, 83, 21, NA, 158, NA, 86, 
                         11, NA, NA, NA, NA, NA, 18, NA, NA, 17, NA, 78, NA, NA, 22, 16, NA, 18, NA, 
                         10, NA, NA, NA, NA, NA, 157, 45, NA, 15, NA, 16, 108, NA, NA, NA, NA, NA, NA, 
                         NA, 102, NA, NA, NA, NA, NA, NA, 18, NA, NA, 19, NA, 12, NA, NA, NA, 156, NA, 
                         67, NA, 58, NA, NA, NA, 36, NA, 18, NA, NA, 11, NA, NA, NA, 8, 34, NA, NA, NA, 
                         NA, NA, 16, NA, NA, 174, 80, NA, 136, 20, NA, NA, NA, NA, NA, NA, NA, 143, NA, 
                         61, NA, NA, NA, 71, NA, NA, NA, NA, NA, NA, 90, NA, 2, NA, NA, NA, NA, NA, NA, 
                         NA, NA, NA, NA, 18, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 94, NA, NA, NA, 
                         151, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 74),.Dim=c(3,56))



#enter times (when individuals are first monitored (days))
S.Y.ENT.S<-structure(.Data=c(6, 4, 5, 2, 5, 0, 0, 5, 5, 4, 5, 5, 5, 0, 3, 5, 5, 6, 3, 0, 5, 4, 5, 2, 6, 5, 
                             3, 5, 5, 3, 2, 4, 5, 3, 4, 6, 4, 1, 2, 5, 2, 5, 2, 4, 4, 0, 4, 4, 3, 3, 4, 5, 
                             1, 3, 5, 6, 3, 3, 6, 2, 5, 5, 2, 4, 2, 2, 1, 5, 1, 5, 3, 4, 1, 3, 5, 3, 0, 5, 
                             4, 5, 6, 4, 6, 4, 4, 5, 4, 6, 3, 5, 5, 1, 2, 5, 5, 5, 4, 2, 2, 5, 3, 6, 1, 4, 
                             0, 4, 5, 4, 4, 3, 5, 6, 4, 6, 4, 5, 4, 2, 2, 5, 5, 4, 5, 4, 5, 3, 2, 4, 0, 2, 
                             4, 5, NA, 4, 4, NA, 3, 0, NA, 5, 5, NA, 3, 4, NA, 3, 5, NA, 6, 4, NA, 3, 0, NA, 
                             4, 6, NA, 5, 4, NA, NA, 0, NA, NA, 2, NA, NA, 4),.Dim=c(3,56))

#censoring times (when they stopped monitoring an individual -- 180 days is the full study period per year)
S.Y.CEN.S<-structure(.Data=c(180, 0, 180, 0, 75, 180, 88, 49, 180, 169, 0, 0, 0, 0, 0, 91, 0, 180, 0, 0, 
                             81, 180, 180, 180, 133, 0, 180, 141, 0, 180, 0, 76, 180, 0, 0, 180, 0, 45, 
                             0, 180, 129, 180, 42, 180, 0, 0, 180, 0, 180, 0, 0, 150, 44, 180, 81, 180, 
                             180, 4, 0, 180, 180, 180, 59, 74, 117, 0, 180, 180, 0, 180, 0, 49, 73, 180, 
                             0, 4, 0, 66, 0, 39, 180, 180, 0, 101, 0, 39, 180, 0, 180, 180, 6, 0, 0, 180, 
                             74, 180, 160, 160, 0, 98, 32, 0, 0, 168, 0, 0, 146, 100, 88, 180, 180, 144,
                             23, 0, 169, 0, 170, 66, 110, 0, 91, 36, 180, 47, 39, 180, 0, 180, 0, 62, 180, 
                             180, NA, 180, 107, NA, 115, 71, NA, 0, 173, NA, 180, 180, NA, 11, 158, NA, 
                             18, 180, NA, 0, 180, NA, 46, 0, NA, 8, 180, NA, NA, 137, NA, NA, 167, NA, NA,
                             0),.Dim=c(3,56))

#number sampled each year
n.obsS<-c(44, 53, 56)

# this is just to visualize the data together a little easier
phi.calf <- data.frame(rbind(S.Y.CEN.S, S.Y.ENT.S, S.Y.S))
rownames(phi.calf) <- c('censor1', 'censor2', 'censor3',
                        'enter1', 'enter2', 'enter3',
                        'exit1', 'exit2', 'exit3')

phi.calf.year1 <- phi.calf[c(1,4,7),]

#####################################################################################################
# winter calf survival data (both sexes)

#exit times
S.Y.W<-structure(.Data=c(53, 50, NA, NA, NA, NA, NA, NA, NA, 10, 33, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                         20, NA, NA, NA, NA, NA, NA, 62, 88, NA, NA, NA, NA, NA, NA, 9, NA, NA, NA, NA, NA, 
                         109, NA, NA, NA, 81, NA, NA, NA, NA, NA, NA, 105, NA, NA, NA, NA, NA, NA, NA, NA, 
                         NA, NA, 92, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 85, NA, NA, NA, 
                         NA),.Dim=c(3,28))

#enter times
S.Y.ENT.W<-structure(.Data=c(0, 0, 0, 5, 0, 0, 4, 0, 0, 4, 0, 0, 3, 0, 0, 4, 0, 0, 3, 0, 0, 4, 0, 0, 4, 
                             0, 0, 4, 0, 0, 4, 0, 0, 4, 0, 0, 4, 2, 0, 4, 2, 0, 5, 2, 0, 0, 2, 0, 0, 3, 
                             0, 0, 3, 0, 0, 3, 0, 0, 3, 0, 0, 3, 0, NA, 3, 0, NA, 3, 0, NA, 3, 0, NA, 3, 
                             NA, NA, 4, NA, NA, 4, NA, NA, 4, NA),.Dim=c(3,28))

#censoring times
S.Y.CEN.W<-structure(.Data=c(0, 0, 185, 185, 148, 185, 185, 54, 185, 0, 0, 184, 185, 122, 185, 185, 111, 
                             183, 185, 185, 185, 0, 104, 185, 185, 165, 185, 185, 0, 0, 48, 2, 184, 185, 
                             107, 184, 0, 185, 185, 185, 185, 185, 0, 185, 185, 16, 0, 185, 74, 185, 181, 
                             4, 185, 0, 185, 185, 185, 34, 185, 183, 36, 185, 185, NA, 0, 183, NA, 175, 
                             185, NA, 9, 182, NA, 185, NA, NA, 185, NA, NA, 0, NA, NA, 119, NA),.Dim=c(3,28))

#number sampled each year
n.obsW<-c(21, 28, 24)

#####################################################################################################
# adult female survival data

#exit times
S.Y.AF<-structure(.Data=c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 75, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                          NA, NA, NA, NA, NA, NA, 161, NA, NA, NA, NA, NA, NA, 187, NA, NA, NA, NA, NA, 
                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                          NA, NA, NA, NA, NA, 320, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                          NA, NA),.Dim=c(3,47))

#enter times
S.Y.ENT.AF<-structure(.Data=c(252, 0, 0, 185, 0, 0, 185, 0, 0, 0, 259, 0, 185, 183, 0, 252, 0, 0, 0, 0, 
                              0, 0, 0, 0, 252, 0, 0, 252, 0, 0, 252, 182, 0, 0, 184, 0, 0, 184, 0, 0, 0, 
                              0, 0, 0, 0, 0, 182, 0, 0, 183, 0, 0, 259, 0, 0, 259, 0, 185, 259, 0, 185, 
                              0, 0, 252, 0, NA, 0, 0, NA, 0, 0, NA, 184, 183, NA, 184, 183, NA, 0, 183, 
                              NA, 252, 259, NA, 252, 259, NA, 252, 0, NA, 0, 0, NA, 0, 0, NA, 0, 259, NA, 
                              0, 259, NA, 0, 259, NA, 0, 0, NA, 0, 184, NA, 184, NA, NA, 184, NA, NA, 0, NA, 
                              NA, 252, NA, NA, 184, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, NA, NA, 0, 
                              NA, NA),.Dim=c(3,47))

#censoring times
S.Y.CEN.AF<-structure(.Data=c(365, 230, 293, 365, 206, 180, 365, 181, 293, 0, 365, 236, 365, 365, 236, 
                              365, 181, 236, 225, 229, 293, 89, 226, 235, 365, 233, 0, 365, 365, 257, 
                              365, 365, 252, 0, 365, 106, 225, 365, 237, 225, 260, 180, 225, 230, 256, 
                              225, 365, 293, 225, 365, 293, 226, 365, 293, 225, 365, 293, 284, 365, 293, 
                              365, 181, 236, 365, 181, NA, 152, 227, NA, 225, 278, NA, 365, 365, NA, 365, 
                              365, NA, 225, 365, NA, 365, 365, NA, 0, 365, NA, 365, 365, NA, 225, 181, NA, 
                              225, 230, NA, 225, 365, NA, 225, 365, NA, 225, 365, NA, 225, 212, NA, 225, 
                              365, NA, 365, NA, NA, 365, NA, NA, 226, NA, NA, 365, NA, NA, 365, NA, NA, 
                              225, NA, NA, 225, NA, NA, 225, NA, NA, 225, NA, NA, 225, NA, NA),.Dim=c(3,47))

#number sampled each year
n.obsAF<-c(47, 37, 21)

#########################################################################################################
# fecundity data

#number pregnant
NumberPregEF<-c(27, 16, 18)

#total sampled
N.TOTALPREGEF<-c(28, 20, 20)

#####################################################################################################
#  aerial count data

#adjusted counts of yearlings (both sexes)
yEF<-c(558, 600, 800, 891)

#adjusted counts of adult females
afEF<-c(2756, 2347, 3190, 2649)

#####################################################################################################
# additional data for model
#number of years included in IPM
nYears=4

#average proportion of yearlings that were female based on ratio of annual male to female calf survival
P=0.584

#initial number of yearlings to start MCMC chains
NinitY<-yEF[1]

#initial number of adult females to start MCMC chains
NinitAF<-afEF[1]

#####################################################################################################
# bundle data for EF (note that counts are transformed to log scale for log-normal distrubtion)
sp.dataEF<-list(nYears = nYears,
                P = P,
                NinitY = NinitY,
                NinitAF = NinitAF,
                log.Y = log(yEF),
                log.AF = log(afEF),
                preg = NumberPregEF,
                afCollars = N.TOTALPREGEF,	
                n.obsAF = n.obsAF, 
                n.obsS = n.obsS, 
                n.obsW = n.obsW,
                yAF = S.Y.AF,
                y.censAF = S.Y.CEN.AF,
                y.entAF = S.Y.ENT.AF,
                yS = S.Y.S,
                y.censS = S.Y.CEN.S,
                y.entS = S.Y.ENT.S,
                yW = S.Y.W,
                y.censW = S.Y.CEN.W,
                y.entW = S.Y.ENT.W)

# End East Fork data				
#####################################################################################################
# Provide initial values for EF for yearling survival (Syguess), shape parameter (shapeS), and
# for indicator variable for dealing with censoring in time-to-event survival models for
# adult females (yAF), summer calf survival (yS), and winter calf survival (yW)

sp.initsEF <- function(){
  Syguess = runif(1, 0, 1)
  shapeguessS=runif(1, 0, 2)
  
  list(
    yAF=with(sp.dataEF, ifelse(is.na(yAF), y.censAF+1, NA)),
    yS=with(sp.dataEF, ifelse(is.na(yS), y.censS+1, NA)),	
    yW=with(sp.dataEF, ifelse(is.na(yW), y.censW+1, NA)),
    shapeS=shapeguessS,	Sy=Syguess)}

#####################################################################################################
# Parmeters to trace in MCMC for EF population. Naming conventions are as follows: 
# yearling survival (Sy), annual calf survival (Sc), summer calf survival (ScS),
# shape parameter for weibull distribution (shapeS), adult female survival (Saf),
# pregnancy probability (Preg; i.e. fecundity), 
# arithmetic mean growth rate (meanGROWTH), geometric mean growth rate (medianGROWTH),
# annual population growth rate (pop.growth; i.e., lambda),number of yearlings (Ny; both sexes),
# number of adult females (Naf), total population size (Ntot), lognormal precison of yearling
# counts #(l.tauY), lognormal precison of adult female counts (l.tauAF), 
# lognormal variance of yearling counts (l.sigma2Y), lognormal variance of adult female counts 
# (l.sigma2AF), lognormal standard deviation of yearling counts (l.sigmaY),
# lognormal standard deviation of adult female counts (l.sigmaAF),
# standard devation of yearling counts on normal scale (sigmaY),
# standard devation of adult female counts on normal scale (sigmaAF)

sp.paramsEF=c("Sy","Sc","ScS","shapeS","ScW","Saf","Preg","meanGROWTH","medianGROWTH","pop.growth",
              "Ny","Naf","Ntot","l.tauY","l.tauAF","l.sigma2Y","l.sigma2AF","l.sigmaY","l.sigmaAF",
              "meanLogNy","meanLogNaf","sigmaY","sigmaAF")

#####################################################################################################
#install and load R2jags packages (note that you will need to install JAGS from website)

#install.packages("R2jags") # run this command if R2jags is not already installed
library(R2jags)

# run base IPM using JAGS (remember to set working directory to folder with JAGS model)

rr.res_EF<-jags(sp.dataEF, sp.initsEF, sp.paramsEF, "JAGS_base_IPM.txt", 
                n.chains=2,  n.iter=200000, n.burnin=150000, n.thin=10)

# end DataS1_EF	   
#####################################################################################################	   

#
##
###
####
#####
####
###
##
#


######################################################################################
# JAGS code for Bitterroot Elk Population IPM
#
######################################################################################
#begin data statement for right censoring indicator
data{
  #Data for exponential adult female survival model-annual
  for(t in 1:(nYears-1)){	
    for(i in 1:(n.obsAF[t])){	
      oneAF[t,i] <- 1  
    }
  }  
  #Data for weibull calf survival model-summer
  for(t in 1:(nYears-1)){	
    for(i in 1:(n.obsS[t])){	
      oneS[t,i] <- 1  
    }
  }  
  #Data for exponential calf survival model-winter
  
  for(t in 1:(nYears-1)){	
    for(i in 1:(n.obsW[t])){	
      oneW[t,i] <- 1  
    }
  }  
}

#start model
model{	
  #Observation model for counts as log normal distribution
  for(t in 1:nYears){
    log.Y[t]~dnorm(log(Ny[t]), l.tauY)
    log.AF[t]~dnorm(log(Naf[t]), l.tauAF)
  }
  #Priors for lognormal precision (tau)
  l.tauY ~ dgamma(0.001, 0.001) 
  l.tauAF ~ dgamma(0.001, 0.001)
  
  #Derive lognormal variance from precision		
  l.sigma2Y<-  1/l.tauY
  l.sigma2AF<- 1/l.tauAF 
  
  #Derive lognormal standard deviation from precision
  l.sigmaY <-  sqrt(l.sigma2Y)
  l.sigmaAF <- sqrt(l.sigma2AF)
  
  #Backtransfrom lognormal standard devation to normal scale
  meanLogNy<-mean(log(Ny))
  meanLogNaf<-mean(log(Naf))
  
  sigmaY <- sqrt(exp(2*meanLogNy+l.sigmaY^2)*(exp(l.sigmaY^2)-1))		
  sigmaAF <- sqrt(exp(2*meanLogNaf+l.sigmaAF^2)*(exp(l.sigmaAF^2)-1))
  
  #priors for shape parameter for summer calf survival
  shapeS ~ dunif(0, 2)
  #informative prior for yearling survival based on Raithel et al. (2007)
  Sy~dbeta(21.92, 2.90)	
  #vague priors for vital rates 
  for(t in 1:(nYears-1)){
    Preg[t]~dbeta(1,1)
    ScS[t]~dbeta(1,1)
    ScW[t]~dbeta(1,1)
    Saf[t]~dbeta(1,1)
    #calculate annual calf survival
    Sc[t]<-ScS[t]*ScW[t]
  }
  ######################################################################################
  # Initial population size (priors) at time 1 as truncated normal distribution
  Ny[1]~dnorm(NinitY, 0.00001)T(0,) 	  #number of yearlings at time 1
  Naf[1]~dnorm(NinitAF, 0.00001)T(0,)  #number of adult females at time 1
  
  ######################################################################################
  # Process model (biological process)
  
  for(t in 2:nYears){
    #Derive number of calves, adults from yearlings, and adults (note that P is proportion of 
    #calves that are female after the first year of life)
    meanNy[t] <- Sc[t-1]*Preg[t-1]*Naf[t-1]
    meanNaf[t] <- (Saf[t-1]*Naf[t-1]) + (Sy*(Ny[t-1]*P))
    
    #process model as Poisson distribution		
    Ny[t] ~ dpois(meanNy[t])
    Naf[t] ~ dpois(meanNaf[t])		
  }		
  for(t in 1:nYears){
    #Derive total population size
    Ntot[t] <- Ny[t] + Naf[t]
    #Derive proportion of population that are yearling 
    propY[t]<-(Ny[t])/(Ntot[t])
    #Derive proportion of population that are adult female 
    propAF[t]<-(Naf[t])/(Ntot[t])
  }
  
  #Derive annual growth rates (with numerical offset)
  for(t in 1:(nYears-1)){
    pop.growth[t]<-((Ntot[t+1] + 0.000001)/(Ntot[t] + 0.000001))
  }
  #Derive arithmetic mean growth rate
  meanGROWTH<- sum(pop.growth)/(nYears-1)
  #Derive geometric mean growth rate
  medianGROWTH<-(prod(pop.growth))^(1/(nYears-1))
  
  ######################################################################################
  #Observation models (aka. Likelihoods for telemetry and pregnancy data)
  #Adult female survival
  for(t in 1:(nYears-1)){
    lambdaAF[t] <- -log(Saf[t])/365
    for(i in 1:(n.obsAF[t])){	
      oneAF[t,i]~dinterval(yAF[t,i], y.censAF[t,i])	
      yAF[t,i]~dexp(lambdaAF[t])T(y.entAF[t,i], )		
    }					
  }
  
  ######################################################################################
  #Pregnancy data from collared adult females 		
  for(t in 1:(nYears-1)){	
    preg[t] ~ dbin(Preg[t], afCollars[t])
  }
  
  ######################################################################################
  #summer calf survival
  
  for(t in 1:(nYears-1)){
    lambdaS[t] <- -log(ScS[t])/pow(185, shapeS)
    for(i in 1:(n.obsS[t])){
      oneS[t,i]~dinterval(yS[t,i], y.censS[t,i])	
      yS[t,i]~dweib(shapeS, lambdaS[t])T(y.entS[t,i], )		
    }			
  }
  
  #winter calf survival
  for(t in 1:(nYears-1)){
    lambdaW[t] <- -log(ScW[t])/180
    for(i in 1:(n.obsW[t])){	
      oneW[t,i]~dinterval(yW[t,i], y.censW[t,i])	
      yW[t,i]~dexp(lambdaW[t])T(y.entW[t,i], )
    }			
  }
}
#End model
######################################################################################