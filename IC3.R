IC3 <- function(Y, Z, b, 
                lambda, gamma,
                Kn, p, D, 
                ku, m, 
                p01, p11, s){
n = nrow(Z)
np = ncol(Z)
eta = Z%*%b
mu = h(eta, p01, p11, s)
Yb = (Y - mu)*h_d(eta, p01, p11, s)/V(mu)
w = h_d(eta, p01, p11, s)^2/V(mu); w = w[,1]
W = sparseMatrix(i = 1:length(w), j = 1:length(w), x = w)
Qma <- as.matrix(bdiag(Penalm(Kn+p, m, lambda[D]),
                       matrix(0, 1+ku, 1+ku)))
for(k in 1:(D-1)){
  Qma <- as.matrix(bdiag(Penalm(Kn+p, m, lambda[D-k]),Qma))
}
Gama <- as.matrix(bdiag(diag((Kn+p)*D)*gamma, matrix(0, 1+ku, 1+ku)))
likelihood = as.numeric(sum(Y*log(mu)+(1-Y)*log(1-mu)))
H = Z%*%solve(t(Z)%*%W%*%Z + Qma + Gama)%*%t(Z)%*%W
traceH = sum(diag(H))
AIC = 2*traceH - 2*likelihood
BIC = log(n)*traceH - 2*likelihood
return(c(likelihood, AIC, BIC))
}