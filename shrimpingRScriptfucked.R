#################################################################
##### Timema Shrimping R Script #################################
#################################################################

# Create Workable Data ----------------------------------------------------

setwd("~/Dropbox/Projects/Shrimping")
shrimpData <- read.table("shrimpData4R.csv", header=T, sep=",")
shrimpData <- shrimpData[,-(27:30)]
shrimpData <- shrimpData[,-18]
shrimpData$unq <- paste(shrimpData$pop2, shrimpData$number)
shrimpData$t1shmp <- ifelse(shrimpData$t1Beh=="shrimp", 1, 0)
shrimpData$t2shmp <- ifelse(shrimpData$t2Beh == "shrimp", 1, 0)
shrimpData$t3shmp <- ifelse(shrimpData$t3Beh == "shrimp", 1, 0)
shrimpData$allShmp1 <- ifelse(shrimpData$t1Beh =="shrimp" |
                             shrimpData$t2Beh == "shrimp" |
                             shrimpData$t3Beh == "shrimp", 1, -1)
shrimpData$twoShmp1 <- ifelse(shrimpData$t1Beh =="shrimp" &
                               shrimpData$t2Beh == "shrimp", 1, -1)
shrimpData$threeShmp <- ifelse(shrimpData$t1Beh =="shrimp" &
                               shrimpData$t2Beh == "shrimp" &
                               shrimpData$t3Beh == "shrimp", 1, 0)
shrimpData$dark1 <- ifelse(shrimpData$color=="red" | shrimpData$color=="grey", 
                          1, -1)
shrimpData$allo1 <- ifelse(shrimpData$population == "hv"|
                          shrimpData$population=="og"|
                          shrimpData$population == "out"| 
                          shrimpData$population == "r12"|
                          shrimpData$population == "log",
                          -1, 1)

expShmp <- data.frame(NULL)

for(i in 1:length(shrimpData$number)) {
  
  expShmp <- rbind(expShmp, shrimpData[i, 1:10], shrimpData[i, 1:10],
                   shrimpData[i,1:10], shrimpData[i, 1:10])
}

expShmp$rep <- rep(1:4, 175)

for(i in 1:length(expShmp$number)) {
if (as.character(shrimpData[expShmp[i,1],expShmp$rep[i]+13]) != "")
  expShmp$beh[i] <- as.character(shrimpData[expShmp[i,1],expShmp$rep[i]+13])
else expShmp$beh[i] <- NA
}
 
for(i in 1:length(expShmp$number)) {
  
    expShmp$beh[i] <- shrimpData[expShmp[i,1],expShmp$rep[i]+10]
}

# Extract Shrimpers -------------------------------------------------------

shmps <- shrimpData[shrimpData$t1Beh=="shrimp" | shrimpData$t2Beh=="shrimp"
                    | shrimpData$t3Beh=="shrimp",]

# Repeatability on two trials
shrimpData1 <- data.frame(shrimpData[-175,])
t1 <- as.vector(shrimpData1$t1Beh)
t2 <- as.vector(shrimpData1$t2Beh)
beh <- c(t1, t2)
beh <- ifelse(beh=="shrimp", 1, 0)
ind <- c(as.vector(shrimpData1$unq), as.vector(shrimpData1$unq))

rpt.mod <- rpt.binomGLMM.multi(beh, ind, "logit", nboot=100, npermut=100)
rpt.mod$P.link # P = 0.01 with 100 perm
rpt.mod$P.org # P = 0.01
rpt.mod$R.link # R = 0.616
rpt.mod$R.org # R = 0.255

###### shrimp GLM ######

shrimpData2 <- shrimpData[-(169:175),]
t1shmpData <-shrimpData2[shrimpData2$t1Beh!="",]
t2shmpData <-shrimpData2[shrimpData2$t2Beh!="",]
t3shmpData <-shrimpData2[shrimpData2$t3Beh!="",]



an1 <- glm(twoShmp ~ dark1, 
           family=binomial, data=t1shmpData, na.action="na.exclude")
summary(an1)
Anova(an1, type="II")

## GLMM ##

for (i in 1:nrow(shrimpData)){
  for (j in 1:ncol(shrimpData)) {
    if (shrimpData[i,j]=="" | is.na(shrimpData[i,j])) {
      shrimpData[i,j] <- NA
    }
  }
}

expShmp <- rbind(shrimpData[,c(1:10, 26, 33, 34)], 
                 shrimpData[,c(1:10, 26, 33, 34)],
                 shrimpData[,c(1:10, 26, 33, 34)])

shmp <- c(shrimpData$t1shmp,shrimpData$t2shmp, shrimpData$t3shmp)
expShmp$shmp <- shmp
expShmp <- na.exclude(expShmp)

library(lme4)
library(MASS)

an1 <- lmer(shmp ~ size + allo1 + (1|unq), family=binomial, data=na.exclude(expShmp))
summary(an1)



an1 <- glm(twoShmp ~ dark1, 
           family=binomial, data=t1shmpData, na.action="na.exclude")
summary(an1)
Anova(an1, type="II")


# plot interaction
library(BradleyTerry2)
library(MASS)
newData <- expand.grid(size=rep(seq(from=min(expShmp$size, na.rm=TRUE), 
                                    to=max(expShmp$size, na.rm=TRUE), 
                                    length.out=100), 2),
                       allo1=c(-1,1))
newData1 <- data.frame(newData, predict(an1, newdata=newData, type="response", se.fit=TRUE))
preds <- predict(an1, newdata=newData, type="response", se.fit=TRUE)

an3 <- MCMCglmm(shmp ~ size, random=~unq, family="categorical", data=expShmp)



