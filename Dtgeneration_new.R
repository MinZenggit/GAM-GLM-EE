##########
# This function generate the simulation dataset
# Ymodel: the outcome models, consider 9 different models, 
# the simulations of Ymodel1,2,5,8 are presented in the main text of our paper
# s01 := P(S=1|Y*=0); s11 := P(S=1|Y*=1)
# p01 := P(Y*=1|Y=0); p11 := P(Y*=1|Y=1)
# n: the source population number.

Dtgeneration_new <- function(Ymodel,
                             emodel,
                             alpha0,
                             p01, p11,
                             s01,
                             s11,
                             n){
  X1 = rnorm(n, 0, 1)
  X2 = runif(n, 0, 1)
  U  = rbinom(n,1,0.5)
  # e been the propensity score of each patient.
  switch(paste0("emodel", emodel),
         emodel1 = {e <- plogis(1 +0.1*X1 -0.1*X2-0.5*U)})
  Tr <- rbinom(n, 1, e)
  #outcome model
  switch(paste0("Ymodel", Ymodel),
         Ymodel1 = {mu = plogis(alpha0 - 2*Tr -1*U - 0.5*X1 + 1*X2)},
         Ymodel2 = {mu = plogis(alpha0 - 2*Tr -1*U - sin(X1*3*pi) +(3*(X2-0.5))^3)},
         Ymodel3 = {mu = plogis(alpha0 - 2*Tr -1*U - exp(2*X1) + (3*(X2-0.5))^3)},
         Ymodel4 = {mu = plogis(alpha0 - 2*Tr -1*U - exp(2*X1) - log(2*X2+0.01))},
         Ymodel5 = {mu = plogis(alpha0 - 2*Tr -1*U - exp(2*X1)- sin(X2*3*pi)*X2)},
         Ymodel6 = {mu = plogis(alpha0 - 2*Tr -1*U - X1^3 - sin(X2*3*pi)*X2)},
         Ymodel7 = {mu = plogis(alpha0 - 2*Tr -1*U - sin(X1*3*pi)*X1 - log(2*X2+0.01))},
         Ymodel8 = {mu = plogis(alpha0 - 2*Tr -1*U - exp(2*X1) + (3*(X2-0.5))^3 + 1*X1*X2)},
         Ymodel9 = {mu = plogis(alpha0 - 2*Tr -1*U  - 0.5*X1 + 1*X2 + 1*X1*X2)})
  Y = rbinom(n,1,mu)
  #generate the misclassified Y
  p_ystar = ifelse(Y, p11, p01)
  Y_star = rbinom(n, 1, p_ystar)
  ps = ifelse(Y_star, s11, s01)
  S = rbinom(n, 1, ps)
  origin=data.frame(Y, X1, X2, U, Tr, Y_star, mu, S, e)
  return(origin)
}
