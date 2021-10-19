
##### Simon's 2 stages #####

install.packages("clinfun")
library(clinfun)

#ph2simon(p0, p1, α, β)
ph2simon(0.2, 0.4, 0.05, 0.1)


####  MCP-MOD procedure #####

install.packages("MCPMod")
library("MCPMod")

##### parameter estimation ####

guesst(d = 0.2, p = 0.9, model = "emax")
guesst(d = c(0.05, 0.2), p = c(0.2, 0.9), model = "logistic")


### model plots ####

models <- list(linear = NULL, emax = c(25),                                
               logistic = c(50, 10.88111), exponential = c(85),            
               betaMod = matrix(c(0.33, 2.31, 1.39, 1.39), byrow=TRUE,nrow=2))

doses <- c(0, 10,25,50,100,150)
plotModels(models, doses, base = 0, maxEff = 0.4, scal = 200,cex=1,lwd=2)


#### the model contrasts and critical value ######


plM <- planMM(models, doses, n = rep(60,6), scal=200, alpha = 0.05)
plM


#### sample size calculation ####

sampSize(models, doses, base = 0, maxEff = 0.4, sigma = 1,             
           upperN = 80, scal = 200, alpha = 0.05)




######### Example 2 #################
# analysing a trial

data(biom)
doses <- c(0, 0.05, 0.2, 0.6, 1)
mods2 <- list(linear = NULL, emax = c(0.05, 0.2), betaMod = c(0.5, 1),
              logistic = matrix(c(0.25, 0.7, 0.09, 0.06), byrow = FALSE,
                                nrow = 2))
plotModels(mods2, doses, base = 0, maxEff = 0.4, scal = 200,lwd=2) 

dfe <- MCPMod(biom, mods2, alpha = 0.05, dePar = 0.05, pVal = TRUE,
              selModel = "maxT", doseEst = "MED2", clinRel = 0.4, off = 1)
summary(dfe)

# plots data with selected model function
plot(dfe, complData = TRUE, cR = TRUE,lwd=2)





