GAMmethod2_wrong <- function(data,
                            Yname,
                            Xname,
                            Uname,
                            Tname,
                            p01, p11, v, s,
                            start ,
                            lambda,
                            gamma ,
                            Kn,
                            p,
                            D,
                            ku,
                            m,
                            iter,
                            epsilon,
                            report){
  n = nrow(data)
  Y_star = data[,Yname]
  Z = as.matrix(data.frame(Tr =  data[,Tname], U = data[,Uname], alpha = rep(1, n)))
  for(Xi in rev(Xname)){
    Zi = bs(data[,Xi], df = Kn+p, degree = p)
    Z = cbind(Zi,Z)
  }
  np = ncol(Z)
  v_star = (p11-p01)*v+p01
  c1 = v_star
  c0 = 1-v_star
  #b estimation
  fit1 = FSAlgo2(Y_star, Z, 
                 start,
                 lambda,gamma,
                 Kn,p,D,ku,m,
                 iter,
                 epsilon,
                 p01, p11, s,
                 report)
  bhat = fit1$bhat
  # tau estimation
  Z1 = Z; Z1[,np-ku] = rep(1, n); eta1 = Z1%*%bhat
  Z0 = Z; Z0[,np-ku] = rep(0, n); eta0 = Z0%*%bhat
  mu11 = mean(plogis(eta1[Y_star==1]))
  mu10 = mean(plogis(eta1[Y_star==0]))
  mu01 = mean(plogis(eta0[Y_star==1]))
  mu00 = mean(plogis(eta0[Y_star==0]))
  thetahats = c(s, bhat, mu11, mu10, mu01, mu00)
  tauhat = c1*mu11+c0*mu10-(c1*mu01+c0*mu00)
  return(tauhat)
}


GLMmethod2_wrong <- function(data,
                            Yname,
                            Xname,
                            Uname,
                            Tname,
                            p01, p11, v, s,
                            start,
                            iter,
                            epsilon,
                            report){
  n = nrow(data)
  Y_star = data[,Yname] 
  Z = as.matrix(data.frame(data[, Xname], Tr =  data[,Tname], data[,Uname], alpha = rep(1, n)))
  np = ncol(Z)
  p01_star = (1-p11)*v/(1-(p11-p01)*v-p01)
  p11_star = p11*v/((p11-p01)*v+p01)
  v_star = (p11-p01)*v+p01
  #s, c, estimation
  c1 = (v-p01_star)/(p11_star-p01_star)
  c0 = (p11_star-v)/(p11_star-p01_star)
  lambda =rep(0, length(Xname));gamma = 0;
  Kn=0;p=1;D=length(Xname);ku = length(Uname)+1;m =0;
  fit1 = FSAlgo2(Y_star, Z, 
                 start,
                 lambda,gamma,
                 Kn,p,D,ku,m,
                 iter,
                 epsilon,
                 p01, p11, s,
                 report)
  bhat = fit1$bhat
  Z1 = Z; Z1[,np-ku] = rep(1, n); eta1 = Z1%*%bhat
  Z0 = Z; Z0[,np-ku] = rep(0, n); eta0 = Z0%*%bhat
  mu11 = mean(plogis(eta1[Y_star==1]))
  mu10 = mean(plogis(eta1[Y_star==0]))
  mu01 = mean(plogis(eta0[Y_star==1]))
  mu00 = mean(plogis(eta0[Y_star==0]))
  thetahats = c(s, bhat, mu11, mu10, mu01, mu00)
  tauhat = c1*mu11+c0*mu10-(c1*mu01+c0*mu00)
  return(tauhat = tauhat)
}