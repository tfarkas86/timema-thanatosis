rm(list=ls())

sd <- read.csv("~/Dropbox/Projects/Shrimping/shmpData2014R.csv")
sd$tmnt <- factor(ifelse(sd$treatment=="once", "control", 
                  ifelse(sd$age=="young", "young", "old")))
sd$tmnt_oVy <- ifelse(sd$tmnt=="old", 1, 0)
sd$tmnt_cVy <- ifelse(sd$tmnt=="control", 1, 0)
#sd <- sd[!sd$treatment=="once",] # only look at Timema with age tmnt
sd$shmp <- ifelse(rowSums(sd[,18:27] > 0), 1, 0)
sd$sex <- ifelse(sd$sex=="M", -.5, .5)
sd$age <- ifelse(sd$age=="young", -.5, .5)
sd <- sd[sd$timID!=13,]

# analysis
contrasts(sd$tmnt) <- cbind(old=c(0,1,0), control=c(1, 0, 0))
an1 <- glm(shmp ~ size + sex + tmnt + phenotype, 
             family=binomial, data=sd) # od = 1.05

an2 <- glmer(shmp ~ size + sex + tmnt + (1|timID), family=binomial, data=sd)
summary(an2)

summary(an1)
Anova(an2, type=3)
lht <- glht(an2, mcp(tmnt="Tukey"))
summary(glht())
summary(lht)

# summary stats
md <- cbind(sd[sd$tmnt=="young",], sd[sd$tmnt=="old",])
t.test(md[,5], md[,36], paired=TRUE)
# plotting data

nd <- data.frame(size=mean(sd$size, na.rm=TRUE), sex=0, age=c(0,1))
nd <- data.frame(nd, predict.glm(an1, newdata=nd, type="response", se.fit=TRUE))
nd$lower <- nd$fit - nd$se.fit; nd$upper <- nd$fit + nd$se.fit

# plotting data with control
an2 <- glm(shmp ~ tmnt_oVy + tmnt_cVy + size + sex + color, family=binomial, data=sd)
summary(an2)
nd <- expand.grid(size=mean(sd$size, na.rm=TRUE), sex=0, 
                  tmnt_oVy=c(0,1), tmnt_cVy=c(0,1))[1:3,]
nd$tmnt <- c("young", "old", "control")
nd <- data.frame(nd, predict.glm(an2, newdata=nd, type="response", se.fit=TRUE))
nd$lower <- nd$fit - nd$se.fit; nd$upper <- nd$fit + nd$se.fit
# plotting

bp <- barplot(height=nd$fit, names.arg=c("young", "old", "control"), 
              ylim=c(0, .7), col="grey", ylab="probability feigning death", 
              las=1)
arrows(x0=bp, y0=nd$lower, x1=bp, y1=nd$upper, angle=90, code=0, lwd=2)

text(bp, rep(.05,3), labels=c("A", "AB", "B"))

# comparison lines
yf <- .05
lines(x=c(bp[1] + 2*(.5/3), bp[1] + 2*(.5/3), bp[2] - xf, bp[2] - xf),
        y = c(nd$fit[1], nd$fit[1] + yf, nd$fit[1] + yf, nd$fit[2]))
lines(x=c(bp[1] + (.5/3), bp[1] + (.5/3), bp[3] + .25, bp[3] + .25),
      y = c(nd$fit[1], nd$fit[1] + 2*yf, nd$fit[1] + 2*yf, nd$fit[3]))
lines(x=c(bp[2] + .25, bp[2] + .25, bp[3] - .25, bp[3] - .25),
      y = c(nd$fit[2], nd$fit[2] + yf, nd$fit[2] + yf, nd$fit[3]))
