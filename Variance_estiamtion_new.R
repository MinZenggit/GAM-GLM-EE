#The total score function
Scoref2_new <- function(Y, Z, 
                        thetas, 
                        lambda, gamma,
                        Kn, p, D, ku, m,
                        p01, p11, v_star){
  n = nrow(Z)
  np = ncol(Z)
  s = thetas[1]
  b = thetas[2:(np+1)]
  mu11 = thetas[np+2]
  mu10 = thetas[np+3]
  mu01 = thetas[np+4]
  mu00 = thetas[np+5]
  #G
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
  G = (t(Z)%*%Yb - Qma%*%b- Gama%*%b)/n
  #mu
  Z1 = Z; Z1[,np-ku] = rep(1, n); eta1 = Z1%*%b
  Z0 = Z; Z0[,np-ku] = rep(0, n); eta0 = Z0%*%b
  Ms = mean(v_star*(1-Y)*s-Y*(1-v_star))
  M1 = mean(Y*(mu11 - plogis(eta1)))
  M2 = mean((1-Y)*(mu10 - plogis(eta1)))
  M3 = mean(Y*(mu01 - plogis(eta0)))
  M4 = mean((1-Y)*(mu00 - plogis(eta0)))
  score = c(Ms, G, M1, M2, M3, M4)
  return(score)
}


#####################################
# the covariance matrix
Getcov2_new <- function(Y, Z, 
                        thetas, 
                        lambda, gamma,
                        Kn, p, D, ku, m,
                        p01, p11, v_star){
  Stheta <- function(thetas){Scoref2_new(Y, Z, 
                                         thetas, 
                                         lambda, gamma,
                                         Kn, p, D, ku, m,
                                         p01, p11, v_star)}
  
  An <-  gradient(Stheta, x = thetas)
  Aninv = solve(An)
  YZ = cbind(Y, Z)
  getphi = function(yz){
    Scoref2_new(yz[1], matrix(yz[-1], nrow = 1),
                thetas,
                rep(0,D), 0,
                Kn, p, D, ku, m,
                p01, p11, v_star) -> phii
    return(phii)
  }
  apply(YZ, 1, getphi) -> phihat
  phihat%*%t(phihat)/nrow(YZ) -> Bn
  cov = Aninv%*%Bn%*%t(Aninv)/nrow(Z)
  return(cov)
}

