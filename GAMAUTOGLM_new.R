GLMmethod2_news <- function(data,
                            Yname,
                            Xname,
                            Uname,
                            Tname,
                            p01, p11, v,
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
  s = mean(Y_star)/v_star/mean(1-Y_star)*(1-v_star)
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
  cov_new = Getcov2_new(Y_star, Z, 
                        thetahats,
                        lambda, gamma,
                        Kn, p, D, ku, m,
                        p01, p11, v_star)
  q_new = c(rep(0, np+1), c1, +c0, -c1, -c0)
  varhat_new = t(q_new)%*%cov_new%*%q_new
  return(list(Y_star = Y_star, Z = Z, tauhat = tauhat, thetahat = bhat, varhat = varhat_new))
}

GAMmethod2_news <- function(data,
                       Yname,
                       Xname,
                       Uname,
                       Tname,
                       p01, p11, v,
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
  s = mean(Y_star)/v_star/mean(1-Y_star)*(1-v_star)
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
  ic = IC3(Y_star, Z, bhat, 
           lambda, gamma,
           Kn, p, D, ku, m, 
           p01, p11, s)
  LK = ic[1]
  AIC = ic[2]
  BIC = ic[3]
  # tau estimation
  Z1 = Z; Z1[,np-ku] = rep(1, n); eta1 = Z1%*%bhat
  Z0 = Z; Z0[,np-ku] = rep(0, n); eta0 = Z0%*%bhat
  mu11 = mean(plogis(eta1[Y_star==1]))
  mu10 = mean(plogis(eta1[Y_star==0]))
  mu01 = mean(plogis(eta0[Y_star==1]))
  mu00 = mean(plogis(eta0[Y_star==0]))
  thetahats = c(s, bhat, mu11, mu10, mu01, mu00)
  tauhat = c1*mu11+c0*mu10-(c1*mu01+c0*mu00)
  cov_new = Getcov2_new(Y_star, Z, 
                        thetahats,
                        lambda, gamma,
                        Kn, p, D, ku, m,
                        p01, p11, v_star)
  q_new = c(rep(0, np+1), c1, +c0, -c1, -c0)
  varhat_new = t(q_new)%*%cov_new%*%q_new
  return(list(Y_star = Y_star, Z = Z, tauhat = tauhat, varhat = varhat_new,
              thetahat = bhat, likelihood = LK, AIC = AIC, BIC = BIC))
}


Gam_auto_news <- function(data,
                     Yname,
                     Xname,
                     Uname,
                     Tname,
                     p01, p11, v,
                     start ,
                     lambdasetting,
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
  v_star = (p11-p01)*v+p01
  c1 = v_star
  c0 = 1-v_star
  np = ncol(Z)
  nlambda = length(lambdasetting)
  AIC = rep(0, nlambda)
  BIC = rep(0, nlambda)
  tauhat = rep(0, nlambda)
  thetahats_matrix = matrix(nrow = nlambda, ncol = np+1+4)
  for(i in 1:nlambda){
    lambda = rep(lambdasetting[i], D)
    GAMmethod2_thetas(data,
                     Yname,
                     Xname,
                     Uname,
                     Tname,
                     p01, p11, v,
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
                     report) %>%try(silent = F) -> fit
    if("try-error" %in% class(fit)){
      AIC[i] = NA
      BIC[i] = NA
      tauhat[i] = NA
      thetahats_matrix[i, ] = NA
    }else{
      AIC[i] = fit$AIC
      BIC[i] = fit$BIC
      tauhat[i] = fit$tauhat
      thetahats_matrix[i, ] = fit$thetahats}
  }
  aic_id = which.min(AIC)
  if(length(aic_id)==1){
    tauhat_aic = tauhat[aic_id]
    thetahats_aic = thetahats_matrix[aic_id, ]
    lambda_aic = rep(lambdasetting[aic_id], D)
    cov_new = Getcov2_new(Y_star, Z, 
                           thetahats_aic,
                           lambda_aic, gamma,
                           Kn, p, D, ku, m,
                           p01, p11, v_star)
    q_new = c(rep(0, np+1), c1, +c0, -c1, -c0)
    varhat_aic = t(q_new)%*%cov_new%*%q_new
  }else{
    tauhat_aic = NA
    varhat_aic = NA
  }
  bic_id = which.min(BIC)
  if(length(bic_id)==1){
    tauhat_bic = tauhat[bic_id]
    thetahats_bic = thetahats_matrix[bic_id, ]
    lambda_bic = rep(lambdasetting[bic_id], D)
    cov_new = Getcov2_new(Y_star, Z, 
                          thetahats_bic,
                          lambda_bic, gamma,
                          Kn, p, D, ku, m,
                          p01, p11, v_star)
    q_new = c(rep(0, np+1), c1, +c0, -c1, -c0)
    varhat_bic = t(q_new)%*%cov_new%*%q_new
  }else{
    tauhat_bic = NA
    varhat_bic = NA
  }
  return(list(tauhat_aic = tauhat_aic, tauhat_bic = tauhat_bic,
              varhat_aic = varhat_aic, varhat_bic = varhat_bic))
}





