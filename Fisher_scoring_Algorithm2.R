#########################
#special link function
h <- function(eta, p01, p11, s){
  return((p01+exp(eta)*p11)*s/(1+p01*(s-1)+exp(eta)*(1+p11*(s-1))))
}
h_d <- function(eta, p01, p11, s){
  exp(eta)*s*(p11-p01)/(1+p01*(s-1)+exp(eta)*(1+p11*(s-1)))^2
}
h_dd <- function(eta, p01, p11, s){
  -exp(eta)*s*(p11-p01)*(exp(eta)*(1+p11*(s-1))-p01*(s-1)-1)/(1+p01*(s-1)+exp(eta)*(1+p11*(s-1)))^3
}
V <- function(mu, p01, p11, s){
  mu*(1-mu)
}
V_d <- function(mu, p01, p11, s){
  1-2*mu
}

########
# The penalty matrix
dfmk <- function(size, order){
  # return a m'th order difference matrix with size of k-m times k
  a = diag(size)
  i = 1
  while(i<=order){
    a = dfm(size-i+1)%*%a
    i = i+1
  }
  return(a)
}
Penalm <- function(size, order, lambda){
  # return the Qm(lambda) metrix
  Delm <- dfmk(size, order)
  return(lambda*t(Delm)%*%Delm)
}
###########
GetGandH2 <- function(Y, Z, b, lambda, gamma,
                      Kn, p, D, 
                      ku, m, 
                      p01, p11, s){
  n = nrow(Z)
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
  G = t(Z)%*%Yb - Qma%*%b- Gama%*%b
  H = -t(Z)%*%W%*%Z-Qma-Gama
  return(list(G = as.vector(G), H))
}
##################
FSAlgo2 <- function(Y,
                   Z,
                   start,
                   lambda,
                   gamma,
                   Kn,
                   p,
                   D,
                   ku,
                   m,
                   iter,
                   epsilon,
                   p01, p11, s,
                   report){
  b0 = start
  bstore = matrix(b0, ncol(Z), 1)
  for(i in 1:iter){
    b_old = as.numeric(bstore[,i])
    GH_old = GetGandH2(Y, Z, b_old, lambda,
                      gamma, Kn, p, D, ku, m,
                      p01, p11, s)
    G_old = GH_old[[1]]
    H_old = GH_old[[2]]
    Gvalue = sum(G_old^2)
    if(report==T){cat(i,"th iteration; 2-norm of G is ",Gvalue,"\n")}
    if(Gvalue < epsilon){
      if(report==T){print("Convergence!")}
      return(list(bhat = b_old, Hessain = H_old))
    }
    b_new = b_old - solve(H_old)%*%G_old
    bstore = cbind(bstore, b_new)
  }
  print("Not convergent!")
}
