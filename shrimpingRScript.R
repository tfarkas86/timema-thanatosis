#################################################################
##### Timema Shrimping R Script #################################
#################################################################

rm(list=ls())

# load and clean data --------------------------------------------

sd <- read.csv("~/Dropbox/Projects/Shrimping/shrimpDataSlim4R.csv", 
               header=T, sep=",")
sd <- sd[sd$pop!="log" & sd$color!="dark",]
sd$c.host <- ifelse(sd$host=="a", -.5, .5)
sd$c.color <- ifelse(sd$color=="striped", -.5, .5)
sd$mal <- ifelse(sd$host=="a", ifelse(sd$color=="green", 1, 0), 
                 ifelse(sd$color=="green", 0, 1))
sd$c.mal <- ifelse(sd$mal==1, .5, -.5)
sd$sex <- ifelse(sd$sex=="male", "male", "female")
sd$c.sex <- ifelse(sd$sex=="male", -.5, .5)
sd$s.size <- scale(sd$size)
sd$s.gflow <- scale(sd$gflow)
sd <- sd[-168,]

# perform statistical analysis -----------------------------------

an1 <- glm(shmp ~ size + sex + host*color*gflow, data=sd, 
             family=binomial)
summary(an1)

an2 <- glm(shmp ~ size + gflow + c.host + c.color + c.sex + c.host * c.color,
           data=sd, family=binomial)
summary(an2)

Anova(an2, type=2)

an3 <- glm(shmp ~ size + gflow + c.host + c.color + c.sex + c.host * c.color,
           data=sd, family=binomial, subset=pop2 != "hvc")
summary(an3)$coefficients[2:7, 1:4]

# mixed model for reviewer
an4 <- glmer(shmp ~ c.host + c.color + s.size  + c.sex + s.gflow + c.host*c.color +
               (1 | pop), family=binomial, data=sd)
summary(an4)

icc <- .5435 / (.5435 + sqrt(var(resid(an4)))) # .35!
phi <- sum(resid(an4, type="pearson")^2)/(df.residual(an4)) # 1.03

write.csv(file="~/Dropbox/Projects/Shrimping/anal_wo_hvc.csv", 
          x=summary(an3)$coefficients)
# create plotting data for size -------------------------------------------

an2 <- glmer(shmp ~ c.host + c.color + c.sex + s.size + s.gflow + (1 | pop), 
           data=sd, family=binomial)
summary(an2)


nd <- expand.grid(s.size=seq(from=min(sd$s.size, na.rm=TRUE), 
                           to=max(sd$s.size, na.rm=TRUE),
                           length.out=100),
                  c.host=0, c.color=0, s.gflow=mean(sd$s.gflow), c.sex=0)

nd <-data.frame(nd, predict(an4, newdata=nd, se.fit=TRUE, type="response"))
nd$upper <- nd$fit + nd$se.fit
nd$lower <- nd$fit - nd$se.fit

# plot prob shrimp against size -----------------------------------
par(mar=c(5,4,4,4))
plot(fit ~ size, data=nd, type="l", ylim=c(-.1, 0.9), xlim=c(10, 24), axes=FALSE,
     xlab="Body length (mm)", ylab="Probability feigning death")
lines(upper ~ size, data=nd, lty=2)
lines(lower ~ size, data=nd, lty=2)

polygon(x=c(nd$size, rev(nd$size)), y=c(nd$lower, rev(nd$upper)), density=200, 
        lwd=.1)

ef <- .1 # width of histogram bars
tr <- .1 # distance below zero for histogram bars

axis(1, pos=0 - tr)
axis(2, at = c(0, .2, .4, .6, .8), labels=c("0.0", "0.2", "0.4", "0.6", "0.8"),
     las=1)
axis(3, labels=FALSE, pos=0.8 + tr) 
axis(4, at=c(-tr, -tr + .05, -tr + .1, -tr + .15, -tr + .2),
     labels=c("0", "5", "10", "15", "20"), pos=24.5, las=1, cex.axis=.8)
mtext(side=4, line=3, text="Sample size", padj=-2, adj=0.05, las=3, cex=.8)
# library(plotrix); axis.break(2, 0.7, style="slash")

zcts <- data.frame(table(sd$size[sd$shmp==0])/100)
colnames(zcts)=c("size","ct")
zcts$size <- 10:23; zcts$shmp <- 0

scts <- data.frame(table(sd$size[sd$shmp==1])/100)
colnames(scts)=c("size","ct")
scts$size <- c(10:17, 19,20); scts$shmp <- .8

for(i in 1:length(zcts$ct)) {
polygon(x=c(zcts$size[i] - ef, zcts$size[i] - ef, zcts$size[i] + ef, 
            zcts$size[i] + ef, zcts$size[i] - ef),
        y=c(zcts$shmp[i] - tr, zcts$ct[i] - tr, zcts$ct[i] - tr, 
            zcts$shmp[i] - tr, zcts$shmp[i] - tr))
}

for(i in 1:length(scts$ct)) {
  polygon(x=c(scts$size[i] - ef, scts$size[i] - ef, scts$size[i] + ef, 
              scts$size[i] + ef, scts$size[i] - ef),
          y=c(scts$shmp[i] + tr, scts$shmp[i]-scts$ct[i] + tr, 
              scts$shmp[i] - scts$ct[i] + tr, scts$shmp[i] + tr, 
              scts$shmp[i] + tr))
}

### Latency

sd <- read.csv("~/Dropbox/Projects/Shrimping/shrimpData_wLat.csv")
sd <- sd[sd$shmp12 == 1,]
sd$lat12avg[sd$lat12avg=="#DIV/0!"] <- NA
sd$dur12avg[sd$dur12avg=="#DIV/0!"] <- NA
sd$lat12avg <- as.numeric(as.character(sd$lat12avg))
sd$dur12avg <- as.numeric(as.character(sd$dur12avg))
hist(sd$lat12avg)
hist(sd$dur12avg)

an1 <- lm(lat12avg ~ pop2, data=sd)
an1 <- lm(dur12avg ~ pop2, data=sd)
Anova(an1)
