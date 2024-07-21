#######################
# Calculate the alpha0 that meets the simulation setted disease prevalence.
# f0 is the disease prevalence
Getalpha0 <- function(n, emodel, Ymodel, f0){
  X1 = rnorm(n, 0, 1)
  X2 = runif(n, 0, 1)
  U  = rbinom(n, 1,0.5)
  # e been the propensity score of each patient.
  switch(paste0("emodel", emodel),
         emodel1 = {e <- plogis(1 +0.1*X1 -0.1*X2-0.5*U)})
  Tr <- rbinom(n, 1, e)
  switch(paste0("Ymodel", Ymodel),
         Ymodel1 = {zb = (-2*Tr -1*U - 0.5*X1 + 1*X2)},
         Ymodel2 = {zb = (-2*Tr -1*U - sin(X1*3*pi)+(3*(X2-0.5))^3)},
         Ymodel3 = {zb = (-2*Tr -1*U - exp(2*X1) + (3*(X2-0.5))^3)},
         Ymodel4 = {zb = (-2*Tr -1*U - exp(2*X1) - log(2*X2+0.01))},
         Ymodel5 = {zb = (-2*Tr -1*U - exp(2*X1)- sin(X2*3*pi)*X2)},
         Ymodel6 = {zb = (-2*Tr -1*U - X1^3 - sin(X2*3*pi)*X2)},
         Ymodel7 = {zb = (-2*Tr -1*U - sin(X1*3*pi)*X1 - log(2*X2+0.01))},
         Ymodel8 = {zb = (-2*Tr -1*U - exp(2*X1) + (3*(X2-0.5))^3 + 1*X1*X2)},
         Ymodel9 = {zb = (-2*Tr -1*U - 0.5*X1 + 1*X2 + 1*X1*X2)})
  fn <- function(alpha0){
    p <- plogis(alpha0 + zb)
    return(mean(p)-f0)
  }
  alpha0 <- uniroot(fn, c(-10, 0))$root
  return(alpha0)
}
a = c()
for(i in 1:100){
  a[i] = a1.XS = Getalpha0(1000000, 1, 1, 0.001)
}

a1.S = Getalpha0(1000000, 1, 1, 0.005)
a1.L = Getalpha0(1000000, 1, 1, 0.01)
a1.M = Getalpha0(1000000, 1, 1, 0.05)
a1.H = Getalpha0(1000000, 1, 1, 0.10)
a2.XS = Getalpha0(1000000, 1, 2, 0.001)
a2.S = Getalpha0(1000000, 1, 2, 0.005)
a2.L = Getalpha0(1000000, 1, 2, 0.01)
a2.M = Getalpha0(1000000, 1, 2, 0.05)
a2.H = Getalpha0(1000000, 1, 2, 0.10)
a3.XS = Getalpha0(1000000, 1, 3, 0.001)
a3.S = Getalpha0(1000000, 1, 3, 0.005)
a3.L = Getalpha0(1000000, 1, 3, 0.01)
a3.M = Getalpha0(1000000, 1, 3, 0.05)
a3.H = Getalpha0(1000000, 1, 3, 0.10)
a4.XS = Getalpha0(1000000, 1, 4, 0.001)
a4.S = Getalpha0(1000000, 1, 4, 0.005)
a4.L = Getalpha0(1000000, 1, 4, 0.01)
a4.M = Getalpha0(1000000, 1, 4, 0.05)
a4.H = Getalpha0(1000000, 1, 4, 0.10)
a5.XS = Getalpha0(1000000, 1, 5, 0.001)
a5.S = Getalpha0(1000000, 1, 5, 0.005)
a5.L = Getalpha0(1000000, 1, 5, 0.01)
a5.M = Getalpha0(1000000, 1, 5, 0.05)
a5.H = Getalpha0(1000000, 1, 5, 0.10)
a6.XS = Getalpha0(1000000, 1, 6, 0.001)
a6.S = Getalpha0(1000000, 1, 6, 0.005)
a6.L = Getalpha0(1000000, 1, 6, 0.01)
a6.M = Getalpha0(1000000, 1, 6, 0.05)
a6.H = Getalpha0(1000000, 1, 6, 0.10)
a7.XS = Getalpha0(1000000, 1, 7, 0.001)
a7.S = Getalpha0(1000000, 1, 7, 0.005)
a7.L = Getalpha0(1000000, 1, 7, 0.01)
a7.M = Getalpha0(1000000, 1, 7, 0.05)
a7.H = Getalpha0(1000000, 1, 7, 0.10)
a8.XS = Getalpha0(1000000, 1, 8, 0.001)
a8.S = Getalpha0(1000000, 1, 8, 0.005)
a8.L = Getalpha0(1000000, 1, 8, 0.01)
a8.M = Getalpha0(1000000, 1, 8, 0.05)
a8.H = Getalpha0(1000000, 1, 8, 0.10)
a9.XS = Getalpha0(1000000, 1, 9, 0.001)
a9.S = Getalpha0(1000000, 1, 9, 0.005)
a9.L = Getalpha0(1000000, 1, 9, 0.01)
a9.M = Getalpha0(1000000, 1, 9, 0.05)
a9.H = Getalpha0(1000000, 1, 9, 0.10)



alpha0list <- c(a1.XS,a1.S,a1.L,a1.M,a1.H,
                a2.XS,a2.S,a2.L,a2.M,a2.H,
                a3.XS,a3.S,a3.L,a3.M,a3.H,
                a4.XS,a4.S,a4.L,a4.M,a4.H,
                a5.XS,a5.S,a5.L,a5.M,a5.H,
                a6.XS,a6.S,a6.L,a6.M,a6.H,
                a7.XS,a7.S,a7.L,a7.M,a7.H,
                a8.XS,a8.S,a8.L,a8.M,a8.H,
                a9.XS,a9.S,a9.L,a9.M,a9.H)
save(alpha0list, file = "True_values_for_simulation/alpha0list.RData")





######################
# Calculate the true ATE values of different simulation settings
Turetau <- function(Ymodel,
                    alpha0,
                    n){
  X1 = rnorm(n)
  X2 = runif(n)
  U  = rbinom(n,1,0.5)
  # e been the propensity score of each patient.
  switch(paste0("Ymodel", Ymodel),
         Ymodel1 = {mu1 = plogis(alpha0 - 2 -1*U - 0.5*X1 + 1*X2);
         mu0 = plogis(alpha0 -1*U - 0.5*X1 + 1*X2)},
         Ymodel2 = {mu1 = plogis(alpha0 - 2 -1*U - sin(X1*3*pi) +(3*(X2-0.5))^3);
         mu0 = plogis(alpha0 -1*U - sin(X1*3*pi) +(3*(X2-0.5))^3)},
         Ymodel3 = {mu1 = plogis(alpha0 - 2 -1*U - exp(2*X1) + (3*(X2-0.5))^3);
         mu0 = plogis(alpha0 -1*U - exp(2*X1) + (3*(X2-0.5))^3)},
         Ymodel4 = {mu1 = plogis(alpha0 - 2 -1*U - exp(2*X1) - log(2*X2+0.01));
         mu0 = plogis(alpha0 -1*U - exp(2*X1) - log(2*X2+0.01))},
         Ymodel5 = {mu1 = plogis(alpha0 - 2 -1*U - exp(2*X1)- sin(X2*3*pi)*X2);
         mu0 = plogis(alpha0 -1*U - exp(2*X1)- sin(X2*3*pi)*X2)},
         Ymodel6 = {mu1 = plogis(alpha0 - 2 -1*U - X1^3 - sin(X2*3*pi)*X2);
         mu0 = plogis(alpha0 -1*U - X1^3 - sin(X2*3*pi)*X2)},
         
         Ymodel7 = {mu1 = plogis(alpha0 - 2 -1*U - sin(X1*3*pi)*X1 - log(2*X2+0.01));
         mu0 = plogis(alpha0 -1*U - sin(X1*3*pi)*X1 - log(2*X2+0.01))},
         
         Ymodel8 = {mu1 = plogis(alpha0 - 2 -1*U  - exp(2*X1) + (3*(X2-0.5))^3 + 1*X1*X2);
         mu0 = plogis(alpha0 -1*U  - exp(2*X1) + (3*(X2-0.5))^3 + 1*X1*X2)},
         
         Ymodel9 = {mu1 = plogis(alpha0 - 2 -1*U  - 0.5*X1 + 1*X2 + 1*X1*X2);
         mu0 = plogis(alpha0 -1*U - 0.5*X1 + 1*X2 + 1*X1*X2)})
  return(mean(mu1-mu0))
}

################
set.seed(25)
n = 1000000
load("True_values_for_simulation/alpha0list.RData")
Ttau1.XS = Turetau(Ymodel = 1,
                   alpha0 = a1.XS,
                   n = n)
Ttau1.S = Turetau(Ymodel = 1,
                  alpha0 = a1.S,
                  n = n)
Ttau1.L = Turetau(Ymodel = 1,
                  alpha0 = a1.L,
                  n = n)
Ttau1.M = Turetau(Ymodel = 1,
                  alpha0 = a1.M,
                  n = n)
Ttau1.H = Turetau(Ymodel = 1,
                  alpha0 = a1.H,
                  n = n)
Ttau2.XS = Turetau(Ymodel = 2,
                   alpha0 = a2.XS,
                   n = n)
Ttau2.S = Turetau(Ymodel = 2,
                  alpha0 = a2.S,
                  n = n)
Ttau2.L = Turetau(Ymodel = 2,
                  alpha0 = a2.L,
                  n = n)
Ttau2.M = Turetau(Ymodel = 2,
                  alpha0 = a2.M,
                  n = n)
Ttau2.H = Turetau(Ymodel = 2,
                  alpha0 = a2.H,
                  n = n)
Ttau3.XS = Turetau(Ymodel = 3,
                   alpha0 = a3.XS,
                   n = n)
Ttau3.S = Turetau(Ymodel = 3,
                  alpha0 = a3.S,
                  n = n)
Ttau3.L = Turetau(Ymodel = 3,
                  alpha0 = a3.L,
                  n = n)
Ttau3.M = Turetau(Ymodel = 3,
                  alpha0 = a3.M,
                  n = n)
Ttau3.H = Turetau(Ymodel = 3,
                  alpha0 = a3.H,
                  n = n)
Ttau4.XS = Turetau(Ymodel = 4,
                   alpha0 = a4.XS,
                   n = n)
Ttau4.S = Turetau(Ymodel = 4,
                  alpha0 = a4.S,
                  n = n)
Ttau4.L = Turetau(Ymodel = 4,
                  alpha0 = a4.L,
                  n = n)
Ttau4.M = Turetau(Ymodel = 4,
                  alpha0 = a4.M,
                  n = n)
Ttau4.H = Turetau(Ymodel = 4,
                  alpha0 = a4.H,
                  n = n)
Ttau5.XS = Turetau(Ymodel = 5,
                   alpha0 = a5.XS,
                   n = n)
Ttau5.S = Turetau(Ymodel = 5,
                  alpha0 = a5.S,
                  n = n)
Ttau5.L = Turetau(Ymodel = 5,
                  alpha0 = a5.L,
                  n = n)
Ttau5.M = Turetau(Ymodel = 5,
                  alpha0 = a5.M,
                  n = n)
Ttau5.H = Turetau(Ymodel = 5,
                  alpha0 = a5.H,
                  n = n)
Ttau6.XS = Turetau(Ymodel = 6,
                   alpha0 = a6.XS,
                   n = n)
Ttau6.S = Turetau(Ymodel = 6,
                  alpha0 = a6.S,
                  n = n)
Ttau6.L = Turetau(Ymodel = 6,
                  alpha0 = a6.L,
                  n = n)
Ttau6.M = Turetau(Ymodel = 6,
                  alpha0 = a6.M,
                  n = n)
Ttau6.H = Turetau(Ymodel = 6,
                  alpha0 = a6.H,
                  n = n)
Ttau7.XS = Turetau(Ymodel = 7,
                   alpha0 = a7.XS,
                   n = n)
Ttau7.S = Turetau(Ymodel = 7,
                  alpha0 = a7.S,
                  n = n)
Ttau7.L = Turetau(Ymodel = 7,
                  alpha0 = a7.L,
                  n = n)
Ttau7.M = Turetau(Ymodel = 7,
                  alpha0 = a7.M,
                  n = n)
Ttau7.H = Turetau(Ymodel = 7,
                  alpha0 = a7.H,
                  n = n)
Ttau8.XS = Turetau(Ymodel = 8,
                   alpha0 = a8.XS,
                   n = n)
Ttau8.S = Turetau(Ymodel = 8,
                  alpha0 = a8.S,
                  n = n)
Ttau8.L = Turetau(Ymodel = 8,
                  alpha0 = a8.L,
                  n = n)
Ttau8.M = Turetau(Ymodel = 8,
                  alpha0 = a8.M,
                  n = n)
Ttau8.H = Turetau(Ymodel = 8,
                  alpha0 = a8.H,
                  n = n)
Ttau9.XS = Turetau(Ymodel = 9,
                   alpha0 = a9.XS,
                   n = n)
Ttau9.S = Turetau(Ymodel = 9,
                  alpha0 = a9.S,
                  n = n)
Ttau9.L = Turetau(Ymodel = 9,
                  alpha0 = a9.L,
                  n = n)
Ttau9.M = Turetau(Ymodel = 9,
                  alpha0 = a9.M,
                  n = n)
Ttau9.H = Turetau(Ymodel = 9,
                  alpha0 = a9.H,
                  n = n)

Ttau = c(Ttau1.XS,Ttau1.S,Ttau1.L,Ttau1.M,Ttau1.H,
         Ttau2.XS,Ttau2.S,Ttau2.L,Ttau2.M,Ttau2.H,
         Ttau3.XS,Ttau3.S,Ttau3.L,Ttau3.M,Ttau3.H,
         Ttau4.XS,Ttau4.S,Ttau4.L,Ttau4.M,Ttau4.H,
         Ttau5.XS,Ttau5.S,Ttau5.L,Ttau5.M,Ttau5.H,
         Ttau6.XS,Ttau6.S,Ttau6.L,Ttau6.M,Ttau6.H,
         Ttau7.XS,Ttau7.S,Ttau7.L,Ttau7.M,Ttau7.H,
         Ttau8.XS,Ttau8.S,Ttau8.L,Ttau8.M,Ttau8.H,
         Ttau9.XS,Ttau9.S,Ttau9.L,Ttau9.M,Ttau9.H)
save(Ttau, file = "True_values_for_simulation/Ttau.RData")

####################
# Truev <- function(Ymodel,
#                   emodel,
#                   alpha0,
#                   design,
#                   n){
#   dt1 <- Dtgeneration(Ymodel,
#                       emodel,
#                       alpha0,
#                       p01 = 0, p11 = 1,
#                       casesize = 1000,
#                       controlsize = 1000,
#                       n)
#   v = mean(dt1$origin$Y)
# }
# n=3000000
# v1.XS = Truev(Ymodel = 1,
#               emodel = 1,
#               alpha0 = a1.XS,
#               n = n)
# v1.S = Truev(Ymodel = 1,
#              emodel = 1,
#              alpha0 = a1.S,
#              n = n)
# v1.L = Truev(Ymodel = 1,
#              emodel = 1,
#              alpha0 = a1.L,
#              n = n)
# v1.M = Truev(Ymodel = 1,
#              emodel = 1,
#              alpha0 = a1.M,
#              n = n)
# v1.H = Truev(Ymodel = 1,
#              emodel = 1,
#              alpha0 = a1.H,
#              n = n)
# 
# v2.XS = Truev(Ymodel = 2,
#               emodel = 1,
#               alpha0 = a2.XS,
#               n = n)
# v2.S = Truev(Ymodel = 2,
#              emodel = 1,
#              alpha0 = a2.S,
#              n = n)
# v2.L = Truev(Ymodel = 2,
#              emodel = 1,
#              alpha0 = a2.L,
#              n = n)
# v2.M = Truev(Ymodel = 2,
#              emodel = 1,
#              alpha0 = a2.M,
#              n = n)
# v2.U = Truev(Ymodel = 2,
#              emodel = 1,
#              alpha0 = a2.H,
#              n = n)
# 
# v3.XS = Truev(Ymodel = 3,
#               emodel = 1,
#               alpha0 = a3.XS,
#               n = n)
# v3.S = Truev(Ymodel = 3,
#              emodel = 1,
#              alpha0 = a3.S,
#              n = n)
# v3.L = Truev(Ymodel = 3,
#              emodel = 1,
#              alpha0 = a3.L,
#              n = n)
# v3.M = Truev(Ymodel = 3,
#              emodel = 1,
#              alpha0 = a3.M,
#              n = n)
# v3.U = Truev(Ymodel = 3,
#              emodel = 1,
#              alpha0 = a3.H,
#              n = n)
# v4.XS = Truev(Ymodel = 4,
#               emodel = 1,
#               alpha0 = a4.XS,
#               n = n)
# v4.S = Truev(Ymodel = 4,
#              emodel = 1,
#              alpha0 = a4.S,
#              n = n)
# v4.L = Truev(Ymodel = 4,
#              emodel = 1,
#              alpha0 = a4.L,
#              n = n)
# v4.M = Truev(Ymodel = 4,
#              emodel = 1,
#              alpha0 = a4.M,
#              n = n)
# v4.U = Truev(Ymodel = 4,
#              emodel = 1,
#              alpha0 = a4.H,
#              n = n)
# v5.XS = Truev(Ymodel = 5,
#               emodel = 1,
#               alpha0 = a5.XS,
#               n = n)
# v5.S = Truev(Ymodel = 5,
#              emodel = 1,
#              alpha0 = a5.S,
#              n = n)
# v5.L = Truev(Ymodel = 5,
#              emodel = 1,
#              alpha0 = a5.L,
#              n = n)
# v5.M = Truev(Ymodel = 5,
#              emodel = 1,
#              alpha0 = a5.M,
#              n = n)
# v5.U = Truev(Ymodel = 5,
#              emodel = 1,
#              alpha0 = a5.H,
#              n = n)
# 
# v6.XS = Truev(Ymodel = 6,
#               emodel = 1,
#               alpha0 = a6.XS,
#               n = n)
# v6.S = Truev(Ymodel = 6,
#              emodel = 1,
#              alpha0 = a6.S,
#              n = n)
# v6.L = Truev(Ymodel = 6,
#              emodel = 1,
#              alpha0 = a6.L,
#              n = n)
# v6.M = Truev(Ymodel = 6,
#              emodel = 1,
#              alpha0 = a6.M,
#              n = n)
# v6.U = Truev(Ymodel = 6,
#              emodel = 1,
#              alpha0 = a6.H,
#              n = n)
# 
# v7.XS = Truev(Ymodel = 7,
#               emodel = 1,
#               alpha0 = a7.XS,
#               n = n)
# v7.S = Truev(Ymodel = 7,
#              emodel = 1,
#              alpha0 = a7.S,
#              n = n)
# v7.L = Truev(Ymodel = 7,
#              emodel = 1,
#              alpha0 = a7.L,
#              n = n)
# v7.M = Truev(Ymodel = 7,
#              emodel = 1,
#              alpha0 = a7.M,
#              n = n)
# v7.U = Truev(Ymodel = 7,
#              emodel = 1,
#              alpha0 = a7.H,
#              n = n)
# 
# v8.XS = Truev(Ymodel = 8,
#               emodel = 1,
#               alpha0 = a8.XS,
#               n = n)
# v8.S = Truev(Ymodel = 8,
#              emodel = 1,
#              alpha0 = a8.S,
#              n = n)
# v8.L = Truev(Ymodel = 8,
#              emodel = 1,
#              alpha0 = a8.L,
#              n = n)
# v8.M = Truev(Ymodel = 8,
#              emodel = 1,
#              alpha0 = a8.M,
#              n = n)
# v8.U = Truev(Ymodel = 8,
#              emodel = 1,
#              alpha0 = a8.H,
#              n = n)
# v9.XS = Truev(Ymodel = 9,
#               emodel = 1,
#               alpha0 = a9.XS,
#               n = n)
# v9.S = Truev(Ymodel = 9,
#              emodel = 1,
#              alpha0 = a9.S,
#              n = n)
# v9.L = Truev(Ymodel = 9,
#              emodel = 1,
#              alpha0 = a9.L,
#              n = n)
# v9.M = Truev(Ymodel = 9,
#              emodel = 1,
#              alpha0 = a9.M,
#              n = n)
# v9.U = Truev(Ymodel = 9,
#              emodel = 1,
#              alpha0 = a9.H,
#              n = n)
# 
