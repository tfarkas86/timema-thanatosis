#################################################################
##### Timema Shrimping R Script #################################
#################################################################

##### Load raw data and manipulate #####
rm(list=ls())
setwd("~/Dropbox/Projects/Shrimping")
shrimpData <- read.table("shrimpData4R.csv", header=T, sep=",")
shrimpData$unq <- paste(shrimpData$pop2, shrimpData$number, sep="")
shrimpData$t1shmp <- ifelse(shrimpData$t1Beh=="shrimp", 1, 
                            ifelse(shrimpData$t1Beh=="", NA, 0))
shrimpData$t2shmp <- ifelse(shrimpData$t2Beh == "shrimp", 1, 
                            ifelse(shrimpData$t2Beh=="", NA, 0))
shrimpData$t3shmp <- ifelse(shrimpData$t3Beh == "shrimp", 1, 
                            ifelse(shrimpData$t3Beh=="", NA, 0))
shrimpData$t4shmp <- ifelse(shrimpData$t4Beh == "shrimp", 1, 
                            ifelse(shrimpData$t4Beh=="", NA, 0))
shrimpData$allShmp1 <- ifelse(shrimpData$t1Beh =="shrimp" |
                             shrimpData$t2Beh == "shrimp" |
                             shrimpData$t3Beh == "shrimp", 1, -1)
shrimpData$twoShmp1 <- ifelse(shrimpData$t1Beh =="shrimp" &
                               shrimpData$t2Beh == "shrimp", 1, -1)
shrimpData$threeShmp <- ifelse(shrimpData$t1Beh =="shrimp" &
                               shrimpData$t2Beh == "shrimp" &
                               shrimpData$t3Beh == "shrimp", 1, 0)
shrimpData$dark1 <- ifelse(shrimpData$color=="red" | 
                             shrimpData$color=="grey", 1, -1)
shrimpData$allo1 <- ifelse(shrimpData$population == "hv"|
                          shrimpData$population=="og"|
                          shrimpData$population == "out"| 
                          shrimpData$population == "r12"|
                          shrimpData$population == "log",
                          -1, 1)
shrimpData$sizeC <- shrimpData$size-mean(shrimpData$size, na.rm=TRUE)
keepCols <- c(3,40,4,31,8,9,10,7,39,32:35, 41)
shpData <- shrimpData[1:174,keepCols]
names(shpData)[4] <- "indv" 

# make long form data farm

shpLong <- rbind(shpData[,c(1:9, 14)], shpData[,c(1:9, 14)], 
                 shpData[,c(1:9, 14)], shpData[,c(1:9, 14)])


for (i in 1:nrow(shpLong)){ # make blanks NA
  for (j in 1:ncol(shpLong)) {
    if (shpLong[i,j]=="" | is.na(shpLong[i,j])) {
      shpLong[i,j] <- NA
    }
  }
}

shrimp1 <- c(shpData$t1shmp,shpData$t2shmp, 
             shpData$t3shmp, shpData$t4shmp)
shpLong$shrimp1 <- shrimp1

###### GLMM analysis ######

## GUIDE TO VARIABLES !
# population: place from which bugs were collected
# allo1: whether population exists in allopatry or sympatry with 
#        populations inhabiting alternate host species; 1 = allopatric
# host: a = Adenostoma, c = Ceanothus
# indv: unique code for an individual, explicit about population
# sex, age, duh
# size: body length
# sizeC: body length, centered around the mean
# color: duh
# dark: green (striped included) vs dark (red and grey); 1 = dark
# shrimp1: whether individual shrimped during trial; 1 = shrimped

## repeatability estimates with rptR library

library(rptR)
repData <- na.omit(shpLong)
# permutations reduced from default (1000) to 10 to speed up
rpt.multi.mod <- rpt.binomGLMM.multi(y=repData$shrimp1, groups=repData$indv, 
                               link="logit", nboot=10, npermut=10)
rpt.add.mod <- rpt.binomGLMM.add(y=repData$shrimp1, groups=repData$indv) 
# the MCMC method (.add) appears not to work with these data, but i DID get it 
# to work before...fucking code. sometimes i think mischeivious elves live in 
# this stupid computer

##GLMM with lmer

library(lme4)
an1 <- lmer(shrimp1 ~ size + allo1 + (1|indv), family=binomial, data=shpLong,
     na.action="na.exclude")
summary(an1)

# with glmmPQL

library(MASS)
an2 <- glmmPQL(shrimp1 ~ sizeC + allo1 + allo1*sizeC, random= ~1|indv, 
               family=binomial, data=shpLong, na.action="na.exclude")
summary(an2)

#### graph PQL resuls: interaction between allo1 and size ####

an3 <- glmmPQL(shrimp1 ~ size + allo1 + allo1*size, random= ~1|indv, 
               family=binomial, data=shpLong, na.action="na.exclude")

newData <- expand.grid(size=seq(from=min(shpLong$size, na.rm=TRUE), 
                                    to=max(shpLong$size, na.rm=TRUE), 
                                    length.out=100),
                       allo1=factor(c(-1,1)))
preds <- predict(an3, newdata=newData, type="response", se.fit=TRUE, level=0)
newData1 <- data.frame(newData, fit=preds)

plot(y=newData1$fit[newData$allo1==-1], x=newData1$size[newData$allo1==-1], type="l", 
     xlab="body length (mm)", ylab="probability of shrimping", ylim=c(0,1))
points(y=newData1$fit[newData$allo1==1], x=newData1$size[newData$allo1==1], type="l", lty=2)
legend("topright", legend=c("allopatric", "parapatric"), bty="n", lty=c(2,1))
