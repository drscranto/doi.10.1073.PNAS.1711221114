##############
## Code to fit the thermal responses of: juvenile mortality, adult mortality, birth rate, and maturation rates for three species:

## We pull data from:

# Lu et al 2010 https://www.jstage.jst.go.jp/article/aez/45/3/45_3_387/_article
# Green Plant bug (Apolygus lucorum) - Temperate species
# data_temp_R.csv

# Amarasekare and Savage 2012 http://www.journals.uchicago.edu/doi/10.1086/663677
# Harlequin bug (Murganta histrionica) - Mediterranean species
# data_med_R.csv

# Dreyer and Baumgartner 1996 http://onlinelibrary.wiley.com/doi/10.1111/j.1570-7458.1996.tb00783.x/abstract
# Pod sucking bug (Clavigralla shadabi) - Tropical species
# data_trop_R.csv

## For each data set we use nls to estimate the parameters of our functions. For some cases where parameters are not estimatable we use a basic least squares function to find parameters that give a local minimum sum of squares within the bounds of biological realism.

## Writes output to a txt file with parameter names and values

## author: Katherine Scranton kscranton@gmail.com
## last updated 9/4/17

##############
rm(list=ls())

# Temperate species:

## ENTER DATA:
data.spp <- read.csv("data_temp_R.csv")

## ESTIMATE PARAMETERS:

# we estimate the birth rate parameters with nls
fecmodel <- nls(brate~bTopt*exp(-((TempK-Toptb)^2)/(2*sb^2)), data=data.spp, start=list(bTopt=2.5, Toptb=298, sb=8))
summary(fecmodel)
write.table(summary(fecmodel)$coefficients[,1], file="temperate_pars.txt",col.names=FALSE)

# we estimate some parameters with nls, but there is insufficient data for colder temperatures to estimate the colder temperature at half maturation rate and its association Arrhenius constant. For that we use a basic least squares function to find a local minimum, enforcing biologically realistic bounds
trm <- 293
mtr <- data.spp[data.spp$TempK==trm,"mat"]
tl <- 273
Al <- -100000
matmodel <- nls(mat ~ mtr*(TempK/trm)*exp(Am*((1/trm) - (1/TempK)))/(1 + exp(Al*((1/tl) - (1/TempK))) + exp(Ah*((1/th) - (1/TempK)))) , data=data.spp, start=list(Am=4000,Ah=35000,th=310))
summary(matmodel)
write.table(summary(matmodel)$coefficients[,1], file="temperate_pars.txt", col.names=FALSE, append=TRUE)
morepars <- c(trm,mtr)
names(morepars) <- c("TRmat","mTR")
write.table(morepars, file="temperate_pars.txt",col.names=FALSE,append=TRUE)

Am <- summary(matmodel)$coefficients[1,1]
Ah <- summary(matmodel)$coefficients[2,1]
th <- summary(matmodel)$coefficients[3,1]

dev.ls <- function(dev.pars){
  Al <- dev.pars[1]
  tl <- dev.pars[2]
  
  diffs <- (data.spp$mat - mtr*(data.spp$TempK/trm)*exp(Am*((1/trm) - (1/data.spp$TempK)))/(1 + exp(Al*((1/tl) - (1/data.spp$TempK))) + exp(Ah*((1/th) - (1/data.spp$TempK)))))^2
  
  sum(diffs,na.rm=TRUE)
}

dev.nls.min <- optim(c(0,291),dev.ls,method="L-BFGS-B", lower=c(-100000,273), upper=c(0,300),control=list(factr=1e-10, pgtol=1e-14))
morepars <- dev.nls.min$par
names(morepars) <- c("ALE","TLE")
write.table(morepars, file="temperate_pars.txt",col.names=FALSE,append=TRUE)

# we estimate the juvenile mortality parameters with nls
tr <- 298
dtr <- data.spp[data.spp$TempK==tr,"jmort"]
jmortmodel <- nls(jmort~dtr*exp(Aj*((1/tr) - (1/TempK))), data=data.spp, start=list(Aj=5000))
summary(jmortmodel)
table.pars <- summary(jmortmodel)$coefficients[,1]
names(table.pars) <- "Adj"
write.table(table.pars, file="temperate_pars.txt",col.names=FALSE, append=TRUE)
morepars <- c(dtr)
names(morepars) <- c("dJTR")
write.table(morepars, file="temperate_pars.txt",col.names=FALSE,append=TRUE)

# we estimate the adult mortality parameters with nls
dtr <- data.spp[data.spp$TempK==tr,"amort"]
amortmodel <- nls(amort~dtr*exp(Aa*((1/tr) - (1/TempK))), data=data.spp, start=list(Aa=5000))
summary(amortmodel)
table.pars <- summary(amortmodel)$coefficients[,1]
names(table.pars) <- "Ada"
write.table(table.pars, file="temperate_pars.txt", col.names=FALSE, append=TRUE)
morepars <- c(tr,dtr)
names(morepars) <- c("TR","dATR")
write.table(morepars, file="temperate_pars.txt",col.names=FALSE,append=TRUE)

#######################
rm(list=ls())

# Mediterranean species:

## ENTER DATA:
data.spp <- read.csv("data_med_R.csv")

## ESTIMATE PARAMETERS:

# we estimate the birth rate parameters with nls
bmodel <- nls(brate~bTopt*exp(-((TempK-Toptb)^2)/(2*sb^2)), data=data.spp, start=list(bTopt=1, Toptb=298, sb=3))
summary(bmodel)
write.table(summary(bmodel)$coefficients[,1], file="mediterranean_pars.txt",col.names=FALSE)

# we estimate some parameters with nls, but there is insufficient data for colder temperatures to estimate the colder temperature at half maturation rate and its association Arrhenius constant. For that we use a basic least squares function to find a local minimum, enforcing biologically realistic bounds
trm <- 297
mtr <- data.spp[data.spp$TempK==trm,"mat"]
tl <- 273
Al <- -100000
matmodel <- nls(mat ~ mtr*(TempK/trm)*exp(Am*((1/trm) - (1/TempK)))/(1 + exp(Al*((1/tl) - (1/TempK))) + exp(Ah*((1/th) - (1/TempK)))) , data=data.spp, start=list(Am=10000,Ah=100000,th=306))
summary(matmodel)
write.table(summary(matmodel)$coefficients[,1], file="mediterranean_pars.txt", col.names=FALSE, append=TRUE)
morepars <- c(trm,mtr)
names(morepars) <- c("TRmat","mTR")
write.table(morepars, file="mediterranean_pars.txt",col.names=FALSE,append=TRUE)

Am <- summary(matmodel)$coefficients[1,1]
Ah <- summary(matmodel)$coefficients[2,1]
th <- summary(matmodel)$coefficients[3,1]

dev.ls <- function(dev.pars){
  Al <- dev.pars[1]
  tl <- dev.pars[2]
  
  diffs <- (data.spp$mat - mtr*(data.spp$TempK/trm)*exp(Am*((1/trm) - (1/data.spp$TempK)))/(1 + exp(Al*((1/tl) - (1/data.spp$TempK))) + exp(Ah*((1/th) - (1/data.spp$TempK)))))^2
  
  sum(diffs,na.rm=TRUE)
}

dev.nls.min <- optim(c(0,295),dev.ls,method="L-BFGS-B", lower=c(-100000,273), upper=c(0,300),control=list(factr=1e-10, pgtol=1e-18))
morepars <- dev.nls.min$par
names(morepars) <- c("ALE","TLE")
write.table(morepars, file="mediterranean_pars.txt",col.names=FALSE,append=TRUE)

# we estimate the juvenile mortality parameters with nls
tr <- 297
dtr <- data.spp[data.spp$TempK==tr,"jmort"]
jmortmodel <- nls(jmort~dtr*exp(Aj*((1/tr) - (1/TempK))), data=data.spp, start=list(Aj=5000))
summary(jmortmodel)
table.pars <- summary(jmortmodel)$coefficients[,1]
names(table.pars) <- "Adj"
write.table(table.pars, file="mediterranean_pars.txt",col.names=FALSE, append=TRUE)
morepars <- c(dtr)
names(morepars) <- c("dJTR")
write.table(morepars, file="mediterranean_pars.txt",col.names=FALSE,append=TRUE)

# we estimate the adult mortality parameters with nls
dtr <- data.spp[data.spp$TempK==tr,"amort"]
amortmodel <- nls(amort~dtr*exp(Aa*((1/tr) - (1/TempK))), data=data.spp, start=list(Aa=5000))
summary(amortmodel)
table.pars <- summary(amortmodel)$coefficients[,1]
names(table.pars) <- "Ada"
write.table(table.pars, file="mediterranean_pars.txt", col.names=FALSE, append=TRUE)
morepars <- c(tr,dtr)
names(morepars) <- c("TR","dATR")
write.table(morepars, file="mediterranean_pars.txt",col.names=FALSE,append=TRUE)

#######################
rm(list=ls())
# Tropical species:

data.spp <- read.csv("data_trop_R.csv")
## ESTIMATE PARAMETERS:

# we estimate the birth rate parameters with nls
bmodel <- nls(brate~bTopt*exp(-((TempK-Toptb)^2)/(2*sb^2)), data=data.spp, start=list(bTopt=10, Toptb=302, sb=1.5))
summary(bmodel)
write.table(summary(bmodel)$coefficients[,1], file="tropical_pars.txt",col.names=FALSE)

# There is insufficient data for extreme temperatures to estimate the colder temperature at half maturation rate and its association Arrhenius constant or any parameters associated with the hot extremes as there is no data for the decay side of the function. So we simply fit a Boltzmann-Arrhenius equation
trm <- 298
mtr <- data.spp[data.spp$TempK==trm,"mat"]
matmodel <- nls(mat ~ mtr*exp(Am*((1/trm) - (1/TempK))), data=data.spp, start=list(Am=5000))
summary(matmodel)
table.pars <- summary(matmodel)$coefficients[,1]
names(table.pars) <- "Am"
write.table(table.pars, file="tropical_pars.txt", col.names=FALSE, append=TRUE)

# we estimate the juvenile mortality parameters with nls
tr <- 298
dtr <- data.spp[data.spp$TempK==tr,"jmort"]
jmortmodel <- nls(jmort~dtr*exp(Aj*((1/tr) - (1/TempK))), data=data.spp, start=list(Aj=5000))
summary(jmortmodel)
table.pars <- summary(jmortmodel)$coefficients[,1]
names(table.pars) <- "Adj"
write.table(table.pars, file="tropical_pars.txt",col.names=FALSE, append=TRUE)
morepars <- c(dtr)
names(morepars) <- c("dJTR")
write.table(morepars, file="tropical_pars.txt",col.names=FALSE,append=TRUE)


# we estimate the adult mortality parameters with nls
dtr <- data.spp[data.spp$TempK==tr,"amort"]
amortmodel <- nls(amort~dtr*exp(Aa*((1/tr) - (1/TempK))), data=data.spp, start=list(Aa=5000))
summary(amortmodel)
table.pars <- summary(amortmodel)$coefficients[,1]
names(table.pars) <- "Ada"
write.table(table.pars, file="tropical_pars.txt", col.names=FALSE, append=TRUE)
morepars <- c(tr,dtr)
names(morepars) <- c("TR","dATR")
write.table(morepars, file="tropical_pars.txt",col.names=FALSE,append=TRUE)


