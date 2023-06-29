rm(list = ls())
library(tidyverse)
library(gmodels)

# 1 BASICS

v1 <- 3
v2 <- c(v1, 1, 2)
v3 <- c(1:5)
v4 <- v2[c(1:2)]
v5 <- c(-1, 2:4, -2)

ma1 <- rbind(v2, 1:3, 2:4)
ma2 <- t(rbind(1:5, 4:8))

v2[1] ; v2[3]  # indexing
sub2 <- ma2[1:2, 1:2]
names(v5)  <- c("a", "b", "c", "d") # name indexing
v5[c("a", "c")]
v5[c(1, 3)]

# labelling matrix
rownames(ma1) <- c("patient1", "patient2", "patient3")
colnames(ma1) <- c("variable1", "variable2", "variable3")

ma3 <- ma1[c("patient1", "patient2")]
ma4 <- ma1[,c("variable1", "variable2")]
# some subset
sub1 <- ma1[c("patient1", "patient2"),c("variable1", "variable3")]

v5[c(T, F, T, F, F)]
v5[v5 < 0] <- 0 # use indexing to set these elements to zero
v5 == 0 # equality
v5 != 0

v5[v5 != 2] <- 1
ma1[ma1[, 1] > 1, ] # IMPORTANT COMMA to ONLY SELECt rows

# LISTS and DF
age <- c(20, 25, 27, 18, 45)
initials <- c("D.T", "H.C", "B.S", "C.B", "M.P")
participants <- data.frame(initials, age, stringsAsFactors = F)
participants[, "age"] ; participants[, "initials"] # indexing in df
participants[participants[, "age"]<25, ] # filter for over 25s
participants$age

# tidyverse funcions
participants %>%
  select(age)  # columns

participants %>%
  filter(age < 25) # boolean search


list1<-list(participants=participants, v2=v2)

# index in lists
list1[[2]]
list1[["v2"]]
list1$v2

# DATASETS
# setwd()

chsi<-read.csv("module0SForStudents/chsiCourse.csv",stringsAsFactors=F) # heatlh data per state
esophCourse <- read.csv("module0SForStudents/esophCourse.csv",stringsAsFactors=F) # histogram data
is.data.frame(chsi) #check
head(chsi)
View(esoph)

chsiCali <- chsi[chsi[, 10] == "California", ]

# STATS

# exploration
hist(chsi$Lung_Cancer) # histogram
boxplot(chsi$Lung_Cancer,chsi$Col_Cancer,chsi$Brst_Cancer) # boxplot
table(chsi$CHSI_State_Name.x) # counts number of occurence of variable in this column
summary(chsi) # stats summary of chsi

# correlations
plot(chsi$Smoker,chsi$Lung_Cancer)
table(esophCourse$esophBn,esophCourse$tobgp) # number of datapoint for these combos
CrossTable(esophCourse$esophBn,esophCourse$tobgp) # some stats with chisqaured

# visualization of contigency table / cross tabulation
ta<-table(esophCourse$esophBn,esophCourse$tobgp) 
mosaicplot(t(ta),col=TRUE)

# ANALYSISING THE DATA

mo<-lm(Lung_Cancer~Smoker,data=chsi)
summary(mo)
plot(Lung_Cancer~Smoker,data=chsi)
abline(mo,col="blue")

# defining functions

multiply<-function(x,y){
  product<-x*y 
  return(product) 
}

multiply(3,4)

circleSurface <- function(r){
  area <- pi * r**2
  return(r)
}

# back to dataanalysis

modSum<-summary(mo)
names(modSum)
modSum$coefficients
modSum$coefficients[,"Pr(>|t|)"] # column of subframe modsum
modSum$coefficients[-1,"Pr(>|t|)"] # removes first row, namely intercept

# p-value functions to extract something from summary table

getPvalues<-function(mo){ #this function takes a model as produced by lm() as an input #make the model summary modSum<-summary(mo) #extract the p values pvalues<-modSum$coefficients[-1,"Pr(>|t|)"] # return the p values as output return(pvalues) }
  #make the model summary
  modSum<-summary(mo)
  pvalues<-modSum$coefficients[-1,"Pr(>|t|)"]
  return(pvalues)
}

getR2<-function(mo){
  modSum<-summary(mo)
  r2<-modSum$r.squared
  return(r2)
}

# call functions
getPvalues(mo)
getR2(mo)


getSSQ<-function(parm){ # input is intercept,slope vector
  # extract slope and intercept
  intercept<-parm[1]
  slope<-parm[2]
  
  yPredicted<-intercept+slope*chsi$Smoker # perfom a linReg fit
  
  lines(chsi$Smoker,yPredicted,col=hsv(runif(1)),lwd=0.25)
  
  squaredDifference<-(yPredicted-chsi$Lung_Cancer)^2 # sum of leats squares per datapoint: measurement of fit-goodness
  
  ssq<-sum(squaredDifference,na.rm=T) # SSS for all datapoints: sum
  
  return(ssq)
}


summary(mo)$coefficients[,"Estimate"]
getSSQ(summary(mo)$coefficients[,"Estimate"]) # input the line of the model 
getSSQ(c(100, -2))  # worse fits
getSSQ(c(25, 10))
getSSQ(c(50, -5))
getSSQ(c(10, 1))

# plot again the datapoints
plot(Lung_Cancer~Smoker,data=chsi)
optimalSSQ<-optim(c(0,0),getSSQ ) # optimal searches in interval 0,0 for optimal function values
opt_par <- optimalSSQ$par
optimalSSQ$value

# Verify that our manual approach with optim and lm() 
# yield the same straight line as the best fit of the data
identical(summary(mo)$coefficients, opt_par) # not identiacal but pretty similar

# some new polynomial fit

getSSQ_quadratic<-function(parm){ # input is intercept,slope vector
  # extract slope and intercept
  a <-parm[1]
  b <-parm[2]
  c <- parm[3]
  
  yPredicted<- a*chsi$Smoker**2 + b*chsi$Smoker + c # perfom a polynomial fit
  
  lines(chsi$Smoker,yPredicted,col=hsv(runif(1)),lwd=0.25)
  
  squaredDifference <- (yPredicted-chsi$Lung_Cancer)^2 # sum of leats squares per datapoint: measurement of fit-goodness
  
  ssq<-sum(squaredDifference,na.rm=T) # SSS for all datapoints: sum
  
  return(ssq)
}

plot(Lung_Cancer~Smoker,data=chsi)
optimalSSQ_quadratic<-optim(c(0,0,0),getSSQ_quadratic) # optimal searches in interval 0,0 for optimal function values
opt_par_q <- optimalSSQ_quadratic$par

# for loops
for(k in 1:5){
  print(k)
}

for(k in 1:5){
  for(kk in 1:5){
    print(c(k,kk))
  }
}

# for loops with if else statements
v6<-c(-1,2,-2,4,5,6,-4); v6positive<-c(); v6negative<-c()

for(k in v6){
  if(k>=0) {
    v6positive<-c(v6positive,k)
  }
  else{
    v6negative<-c(v6negative,k) 
  }
}

# same loop but with booleans

v6<-c(-1,2,-2,4,5,6,-4); v6index<-c()

for(k in v6){
  if(k>=0) {
    v6index<-c(v6index,T)
  }
  else{
    v6index<-c(v6index,F) 
  }
}

v6pos <- v6[v6index]
v6neg <- v6[!v6index]

# systematically asses all ssess the association of potential exposures with the incidence of
# health-related outcomes in the CHSI data:

colnames(chsi)

R2matrix<-matrix(NA,nrow=5,ncol=4) # matrix to be filled with R2
Pmatrix<-matrix(NA,nrow=5,ncol=4) # with P

rownames(R2matrix)<-colnames(chsi)[1:5] # sqaure matrix displaying associations of columns
colnames(R2matrix)<-colnames(chsi)[6:9]

for(k in 1:5){
  for(j in 6:9){
    mo<-lm(paste((colnames(chsi)[k]),"~",(colnames(chsi)[j])),data=chsi)
    R2matrix[k,j-5]<-getR2(mo)
    Pmatrix[k,j-5]<-log10(getPvalues(mo))
  }
}

heatmap(R2matrix,col = gray.colors(start=0.9, end=0.2,256),scale = "none",Rowv =NA,Colv=NA) # generates heatmap of this matrix
heatmap(Pmatrix,col = gray.colors(start=0.9, end=0.2,256),scale = "none",Rowv =NA,Colv=NA) # generates heatmap of this matrix

print("Causility and Correlaion of features, we see what correaltes")

# logistic regression with glm in esophCourse


# POPULATION DYNAMICS

logisticMapDynamics<-function(xInitial,growthRate,nGenerations,nBurnIn){
  dynamics<-rep(NA,nGenerations)
  dynamics[1]<-xInitial
  for(k in 1:(nGenerations-1)){
    dynamics[k+1]<-growthRate*dynamics[k]*(1-dynamics[k])
  }
  dynamics[nBurnIn:nGenerations]
}

plot(logisticMapDynamics(0.01,1.5,100,10) ,type="l",ylim=c(0,1))
lines(logisticMapDynamics(0.01,2.5,100,10),col=3)

res<-c() 
for(r in seq(1.3,3.9,by=0.005)){
  res<-cbind(res,rbind(r, logisticMapDynamics(0.01,r,500,100))) 
}

plot(res[1,],res[2,],pch=".") # emerging chaos

logisticMapDynamics2<-function(xInitial,growthRate,nGenerations,nBurnIn){
  
  dynamics<-rep(NA,nGenerations)
  dynamics[1]<-xInitial
  
  for(k in 1:(nGenerations-1)){
    dynamics[k+1]<-dynamics[k]*exp(growthRate * (1-dynamics[k])) # advantage is that pop doesn't become negative
  }
  dynamics[nBurnIn:nGenerations]
}

res1<-c() 
for(r in seq(1.3,3.9,by=0.005)){
  res1<-cbind(res1,rbind(r, logisticMapDynamics(0.01,r,500,100))) 
}

plot(res1[1,],res1[2,],pch=".")
