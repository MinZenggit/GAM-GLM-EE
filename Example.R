# An illustration of methods

library(Matrix)
library(extRC)
library(splines)
library(dplyr)
library(rootSolve)
rm(list = ls())
load("True_values_for_simulation/alpha0list.RData") # True alpha0 values
load("True_values_for_simulation/Ttau.RData") # True ate values
source("Dtgeneration_new.R")                  # data generation function
source("GAMAUTOGLM_new.R")                    # the main GAM,GLM methods function
source("Variance_estiamtion_new.R")           # variance estimation function for GAM/GLM methods 
source("Fisher_scoring_Algorithm2.R")         # optimization algorithm for GAM/GLM methods 
source("IC3.R")                               # Information criterion calculation function for GAM method
source("simulationsettings.R")                # simulation settings

#############################
# Generate an example data set
Ymodel = 1 # Ymodel1 = {mu = plogis(alpha0 - 2*Tr -1*U - 0.5*X1 + 1*X2)}
emodel = 1 # emodel1 = {e <- plogis(1 +0.1*X1 -0.1*X2-0.5*U)}
p01 = 0; p11 = 0.8;# p01 := P(Y*=1|Y=0); p11 := P(Y*=1|Y=1)
v = 0.01 # disease prevalence
alpha0 = -6.306 # calculate to make v = 0.001
casesize =2000; controlsize =2000; n =2000000;
# casesize/controlsize: the size of selected case and
# control samples from source population
Trau = -0.0206245315 # the true ate value
gets(n, casesize, controlsize, p01, p11, v = v)-> a 
# According the case-control size to calculate 
# the selection probabilities for data generation.
# s01 := P(S=1|Y*=0); s11 := P(S=1|Y*=1)
s01 = a[1]; s11 = a[2]

dt1 <- Dtgeneration_new(Ymodel,
                        emodel,
                        alpha0,
                        p01, p11,
                        s01,
                        s11,
                        n) 
dt1 %>% filter(S==1)%>% select(c("X1", "X2", "U", "Tr", "Y_star")) -> data1
# data1 is the case-control dataset

###########################
# Apply the GAM/GEM-EE methods.
Kn = 10 # knots of spline
p = 3; # degree of spline
D = 2; # number of continous covariates
ku = 2; # number of category covariates
m = 2 # order of difference matrix
lambda = c(5, 5); gamma = 0.01
# Tuning parameters
iter = 200; epsilon = 0.000001
# ends points

GAMmethod2_news(data1,
                Yname = "Y_star", # Outcome variable
                Xname = c("X1", "X2"), # continous covariates
                Uname = "U", # category covariates
                Tname = "Tr", # treatment variable
                p01, p11, v,
                start = c(rep(0, D*(Kn+p)+ku), -3), # initial values 
                lambda,
                gamma ,
                Kn,
                p,
                D,
                ku,
                m,
                iter,
                epsilon,
                report) %>% try(silent = F) -> fit1
if("try-error" %in% class(fit1)){
  tauhat_gam = NA; varhat_gam = NA
}else{
  tauhat_gam = fit1$tau;
  varhat_gam = as.numeric(fit1$varhat)}
GLMmethod2_news(data1,
                Yname = "Y_star", # Outcome variable
                Xname = c("X1", "X2"), # continous covariates
                Uname = "U", # category covariates
                Tname = "Tr", # treatment variable
                p01, p11, v,
                start = c(0, 0, 0, 0, -3),
                iter,
                epsilon,
                report) %>% try(silent = F) -> fit2
if("try-error" %in% class(fit2)){
  tauhat_glm  = NA; varhat_glm = NA
}else{
  tauhat_glm = fit2$tau;
  varhat_glm = as.numeric(fit2$varhat)}

c(tauhat_gam, varhat_gam, tauhat_glm, varhat_glm)
# True ate: Trau = -0.0206245315 

