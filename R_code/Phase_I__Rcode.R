
##### CRM, mCRM #####
install.packages("bcrm")
library("bcrm")


#patient = 1:42
dose.label <- c(5, 10, 15, 25, 40, 50, 60)
p.tox0 <- c(0.05, 0.10, 0.20, 0.30, 0.35, 0.40, 0.45)

#### CRM ###
trial.crm <- bcrm(stop = list(nmax = 12), p.tox0 = p.tox0, dose = dose.label,
                     ff = "ht", prior.alpha = c(1, 1, 1), target.tox = 0.33,
                     constrain = FALSE, plot=TRUE, cohort = 1,pointest = "plugin",
                     method="exact") 
                               
plot(trial.crm,trajectory = TRUE)
print(trial.crm)


#### mCRM ###
trial.mcrm <- bcrm(stop = list(nmax = 24), p.tox0 = p.tox0, dose = dose.label,
                     ff = "ht", prior.alpha = c(1, 1, 1), target.tox = 0.33,
                     constrain = TRUE,plot=TRUE,cohort = 3,pointest = "plugin", 
                     method="exact",start=1)          


plot(trial.mcrm,trajectory = TRUE)
print(trial.mcrm)



##### TPI, mTPI, mTPI-2 ####


install.packages("remotes")
remotes::install_github("anirban0451/tpidesigns")

library(tpidesigns)


#Plotting of Posterior Distribution and UPM for TPI
n0 <- 3  ######
x0 <- 0  ######

UPMTPI <- function(n0,x0) {
  n0 <- n0
  x0 <- x0
  alp <- 0.005+x0
  beta <- 0.005+n0-x0
  sigma <- sqrt((alp*beta)/(((alp+beta)^2)*(alp+beta+1)))
  pt <- 0.3
  e10 = 1.5*sigma
  e20 = 1*sigma
  if (e10 > pt) {e10 <- pt-0.005}
  if (e20 > (1-pt)) {e20 <- 1-pt}
  
  
  p1 <- pbeta(pt-e10,alp,beta)
  p2 <- pbeta(pt+e20,alp,beta)
  pe <- round(p1,4)
  ps <- round(p2-p1,4)
  pd <- round(1-p2,4)
  pp <- round(c(pe,ps,pd),4)
  #pp
  qq <- c(paste("p(E)=",pe),paste("p(S)=",ps),paste("p(D)=",pd))
  print(qq)
  #print(qq,pt-e10,pt+e20)
  
  
  upmplot(x = x0, n = n0, pt = 0.3, design = "mtpi",e1 = e10, e2 = e20,  
          w=0,a1 = 0.005, b1 = 0.005, a2 = alp, b2 = beta)
  
}

UPMTPI(3,0)
UPMTPI(3,1)
UPMTPI(3,2)
UPMTPI(3,3)



#Plotting of Posterior Distribution and UPM for TPI
x0 <- 4  ######
n0 <- 9  ######

alp <- 0.005+x0
beta <- 0.005+n0-x0
sigma <- sqrt((alp*beta)/(((alp+beta)^2)*(alp+beta+1)))
pt <- 0.3
e10 = 1.5*sigma
e20 = 1*sigma
if (e10 > pt) {e10 <- pt-0.005}
if (e20 > (1-pt)) {e20 <- 1-pt}


p1 <- pbeta(pt-e10,alp,beta)
p2 <- pbeta(pt+e20,alp,beta)
pe <- round(p1,4)
ps <- round(p2-p1,4)
pd <- round(1-p2,4)
pp <- round(c(pe,ps,pd),4)
pp
qq <- c(paste("p(E)=",pe),paste("p(S)=",ps),paste("p(D)=",pd))
list(print(qq),pt-e10,pt+e20)

upmplot(x = x0, n = n0, pt = 0.3, design = "mtpi",e1 = e10, e2 = e20,  
        w=0,a1 = 0.005, b1 = 0.005, a2 = alp, b2 = beta)

pt-e10
pt+e20



#Plotting of Posterior Distribution and UPM for mTPI
x0 <- 4  ######
n0 <- 9  ######
alp <- 0.005+x0
beta <- 0.005+n0-x0

pt <- 0.3
e10 = 0.05
e20 = 0.05

p1 <- pbeta(pt-e10,alp,beta)
p2 <- pbeta(pt+e20,alp,beta)
pe <- round(p1/(pt-e10),4)
ps <- round((p2-p1)/(e10+e20),4)
pd <- round((1-p2)/(1-pt-e20),4)
qq <- c(paste("UPM(E)=",pe),paste("UPM(S)=",ps),paste("UPM(D)=",pd))
print(qq)

upmplot(x = x0, n = n0, pt = 0.3, design = "mtpi",e1 = e10, e2 = e20,  
        w=0,a1 = 0.005, b1 = 0.005, a2 = alp, b2 = beta)
pt-e10
pt+e20



#Plotting of Posterior Distribution and UPM for mTPI-2
x0 <- 4  ######
n0 <- 9  ######
alp <- 0.005+x0
beta <- 0.005+n0-x0

pt <- 0.3
e10 = 0.05
e20 = 0.05

p1 <- pbeta(pt-e10,alp,beta)
p2 <- pbeta(pt+e20,alp,beta)
pe <- round(p1/(pt-e10),4)
ps <- round((p2-p1)/(e10+e20),4)
pd <- round((1-p2)/(1-pt-e20),4)
qq <- c(paste("UPM(E)=",pe),paste("UPM(S)=",ps),paste("UPM(D)=",pd))
#print(qq)

#upmplot(x = x0, n = n0, pt = 0.3, design = "mmtpi",e1 = e10, e2 = e20,  
#        w=0,a1 = 0.005, b1 = 0.005, a2 = alp, b2 = beta)

upmplot(x = x0, n = n0, pt = 0.3, design = "mmtpi",e1 = e10, e2 = e20,  
        w=0,a1 = 0.005, b1 = 0.005, a2 = alp, b2 = beta)
pt-e10
pt+e20


####### BOIN #####

install.packages("BOIN")
library("BOIN")


#### Obtain dose escalation and de-escalation boundaries ####
bound <- get.boundary(target = 0.3, ncohort = 15, cohortsize = 2)
summary(bound)

plot(bound)



### select MTD ####

n <- c(1,1,8,17)
y <- c(0,0,1,5)
selmtd <- select.mtd(target=0.3, npts=n, ntox=y)
summary(selmtd)
plot(selmtd)






