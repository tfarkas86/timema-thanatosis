rm(list=ls())

setwd("~/Dropbox/Projects/Shrimping")

fdw <- read.table("shrimpData4R.csv", header=TRUE, sep=",") # read all data
fdw <- fdw[order(fdw$no),] # sort by individual
drop <- names(fdw)[c(2,11:13, 18:30)] # list of columns to remove
drop2 <- names(fdw)[c(2,11:17, 18, 19, 21, 23, 25, 27:30)] # list of columns to remove

fdl <- reshape(data=fdw, varying=list(c("t1Beh", "t2Beh", "t3Beh", "t4Beh")), 
                direction="long", drop=drop, timevar="trial") # long form
fdl2 <- reshape(data=fdw, varying=list(c("t1Dur", "t2Dur", "t3Dur", "t4Dur")), 
               direction="long", drop=drop2, timevar="trial") # long form
fdl2$t1Dur <- ifelse(is.na(fdl2$t1Dur), 0, fdl2$t1Dur )

fdl$shmp <- ifelse(fdl$t1Beh=="shrimp", 1, ifelse(fdl$t1Beh=="" ,NA , 0)) # code shrimping as 1
fdl$allo <- factor(ifelse(fdl$population=="pr"|fdl$population=="s", "A", "P")) # code allopatric vs. parapatric
fdl$dark <- factor(ifelse(fdl$color=="red"|fdl$color=="grey", "dark", "green"))
fdl$gflow <- if(fdl$population=="s" | fdl$population=="pr") 100
else if (fdl$pop2=="ogc") .99
else if (fdl$pop2=="oga"), .1
else if (fdl$pop2=="r12c"), 
                                                                            ))

fdl2$shmp <- ifelse(fdl$t1Beh=="shrimp", 1, ifelse(fdl$t1Beh=="" ,NA , 0)) # code shrimping as 1
fdl2$allo <- factor(ifelse(fdl$population=="pr"|fdl$population=="s", "A", "P")) # code allopatric vs. parapatric
fdl2$dark <- factor(ifelse(fdl$color=="red"|fdl$color=="grey", "dark", "green"))


# binomial models

mod1 <- lmer(shmp ~ size + sex + age + dark + host + allo + (1|id), 
             data=fdl, family=binomial); summary(mod1)
mod2 <- lmer(shmp ~ size + allo + (1|id), data=fdl, subset=fdl$trial<5, family=binomial, 
             na.rm=TRUE); summary(mod2)
