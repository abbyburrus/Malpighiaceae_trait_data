################Processing tree to fit what coding was available in literature######################
library(ape)
library(geiger)
library(phytools)
library(ggtree)
install.packages("lmtest")
library(ggplot2)
library(phangorn)
library(zoo)
library(lmtest)
read.tree("2015tree.tre")->malpigh
#this is just a slightly-trimmed version of the 2015 Malpighiaceae phylogeny

statestable <- read.table("shortenedtipssingle.csv", header = TRUE, sep = ",", colClasses='character')
#shortenedtips.csv is just a .csv version of the matrix trimmed to the taxa that are present in the 2015 Malpighiaceae phylogeny 

#read in your traits, one at a time (or, at least, each separated into their own set of values/character string)
states1 <- as.vector(statestable$sepalgland_phytools)
names(states1)<-statestable$xx

states1 <- as.vector(statestable$sepalgland_phytools)
states2 <- as.vector(statestable$lateral_organ) 
states3 <- as.vector(statestable$floralnec)

names(states1)<-statestable$xx
names(states2)<-statestable$xx
names(states3)<-statestable$xx
states4 <- as.vector(statestable$petal)
states5 <- as.vector(statestable$anther) 

names(states4)<-statestable$xx
names(states5)<-statestable$xx

states6 <- as.vector(statestable$lam_blade_phytools)
states7 <- as.vector(statestable$petiole_phytools)
names(states6)<-statestable$xx
names(states7)<-statestable$xx

states8 <- as.vector(statestable$NWOW)
names(states8)<-statestable$xx

states9 <- as.vector(statestable$calyx)
names(states9)<-statestable$xx


####################Simmap demonstration with sepal glands:
sepaltest<-make.simmap(malpigh, states1, nsim=100) #default model is ARD (all rates different) - if this is no longer the case when you work with phytools, include the parameter model="ARD"
sumsepaltest<-summary(sepaltest)

binarysepaltest<-make.simmap(malpigh, states9, nsim=100) #this is the same, but for the binary sepal condition of glandular or eglandular
sepalbin<-summary(binarysepaltest)

####################Plotting demonstration with this example:
cols<-setNames(c("black","red","orange","yellow","green"),levels(c("0","1","2","3","4")))
plotTree(malpigh,fsize=0.5,ftype="i",lwd=2, offset=0.8)
nodelabels(node=1:malpigh$Nnode+Ntip(malpigh),
           pie=sumsepal$ace,piecol=cols,cex=0.08)
tiplabels(pie=sumsepal$tips,piecol=cols,cex=0.025)

#If working with a binary character, specify binary colors to represent the states:
cols2<-setNames(c("black","pink"),levels(c("0","1")))
plotTree(malpigh,fsize=0.5,ftype="i",lwd=2, offset=0.8)
nodelabels(node=1:malpigh$Nnode+Ntip(malpigh),
           pie=sepalbin$ace,piecol=cols2,cex=0.08)
tiplabels(pie=sepalbin$tips,piecol=cols2,cex=0.025)
#save as image (.png) with image dimensions 2000x2000 pixels- bigger or smaller depending on your own phylogeny


####################Incorporating a multistate set of models to see if various traits change more/less in the OW vs. NW or not

###regime map of a single stochastic mapping, incorporating tip data/simulated history and tip data of OW vs. NW 
regionmap<-make.simmap(malpigh, states8, model="ARD")
regionmap2<-make.simmap(malpigh, states8, model="ARD") 
regionmap3<-make.simmap(malpigh, states8, model="ARD")
#The stochasticity of the single SIMMAP generation doesn't seem to change the outcome- the level of significance seems about the same across the board of stochastic mappings.
fit.multi<-fitmultiMk(regionmap2,states1)
#your trait needs to be a distinct vector of its own- but I already gave it that by how I imported the data
print(fit.multi,digits=2)

#buiding single-state model for contrast, letting us test whether the multistate model or single-state model better explains the data:
fit.single<-fitMk(regionmap2,states1,model="ER") 
print(fit.single,digits=2)
#Both have very similar log-likelihoods

lrtest(fit.single,fit.multi)
#This gives results. Example output look something like this:
#sepal glands according to region:
#Model 1: fit.single
#Model 2: fit.multi
#Df  LogLik Df Chisq Pr(>Chisq)
#1   1 -421.04                    
#2   2 -417.94  1 2.195     0.0385
#The generally accepted 0.05 p-value is applied here. 


####################Incorporating Pagel's 1994 model of correlation to test for character nonindependence
fitPagel(malpigh,states1,states2)->callam
fit.bcsepdep<-fitPagel(malpigh,states1,states2,dep.var="x")
fit.bclatdep<-fitPagel(malpigh,states1,states2,dep.var="y")

elailamaic<-setNames(c(callam$independent.AIC,
                       callam$dependent.AIC,
                       fit.bcsepdep$dependent.AIC,
                       fit.bclatdep$dependent.AIC),
                     c("independent","dependent x&y",
                       "dependent sepal glands","dependent leaf glands"))
elailamaic
aic.w(elailamaic)
plot(callam,lwd.by.rate=TRUE)

#> elailamaic
#independent          dependent x&y dependent sepal glands  dependent leaf glands 
#640.7566               604.2908               622.4070               610.6549 
#> aic.w(elailamaic)
#independent          dependent x&y dependent sepal glands  dependent leaf glands 
#0.00000001             0.96004727             0.00011179             0.03984093 
#Interdependence (Non-independence)
#In this example, the interdependent model seems to explain the data best.

