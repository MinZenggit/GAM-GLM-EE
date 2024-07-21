simu_single_new <- function(Ymodel,
                            Yname,
                            Xname,
                            Uname,
                            Tname,
                            emodel,
                            alpha0,
                            p01, p11,
                            s01, s11,
                            n,
                            v,
                            startgam,
                            startglm,
                            lambda,
                            gamma,
                            Kn,
                            p,
                            D,
                            ku ,
                            m,
                            iter,
                            epsilon,
                            report){
  dt1 <- Dtgeneration_new(Ymodel,
                          emodel,
                          alpha0,
                          p01, p11,
                          s01,
                          s11,
                          n) 
  dt1 %>% filter(S==1) -> data1
  rm(dt1)
  GAMmethod2_news(data1,
             Yname,
             Xname,
             Uname,
             Tname,
             p01, p11, v,
             start = startgam,
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
                 Yname,
                 Xname,
                 Uname,
                 Tname,
                 p01, p11, v,
                 start = startglm,
                 iter,
                 epsilon,
                 report) %>% try(silent = F) -> fit2
  if("try-error" %in% class(fit2)){
    tauhat_glm  = NA; varhat_glm = NA
  }else{
    tauhat_glm = fit2$tau;
    varhat_glm = as.numeric(fit2$varhat)}
  return(c(tauhat_gam, varhat_gam, tauhat_glm, varhat_glm))
}


simu_auto_new<- function(Ymodel,
                        Yname,
                        Xname,
                        Uname,
                        Tname,
                        emodel,
                        alpha0,
                        p01, p11,
                        s01, s11,
                        n,
                        v,
                        startgam,
                        startglm,
                        lambdasetting,
                        gamma,
                        Kn,
                        p,
                        D,
                        ku ,
                        m,
                        iter,
                        epsilon,
                        report){
  dt1 <- Dtgeneration_new(Ymodel,
                          emodel,
                          alpha0,
                          p01, p11,
                          s01,
                          s11,
                          n) 
  dt1 %>% filter(S==1) -> data1
  Gam_auto_news(data1,
           Yname,
           Xname,
           Uname,
           Tname,
           p01, p11, v,
           start = startgam,
           lambdasetting,
           gamma ,
           Kn,
           p,
           D,
           ku,
           m,
           iter,
           epsilon,
           report) %>% try(silent = F)-> fit1
  if("try-error" %in% class(fit1)){
    tauhat_gam_aic = NA
    tauhat_gam_bic = NA
    varhat_gam_aic = NA
    varhat_gam_bic = NA
  }else{
    tauhat_gam_aic = fit1$tauhat_aic
    tauhat_gam_bic = fit1$tauhat_bic
    varhat_gam_aic = fit1$varhat_aic
    varhat_gam_bic = fit1$varhat_bic}
  GLMmethod2_news(data1,
                  Yname,
                  Xname,
                  Uname,
                  Tname,
                  p01, p11, v,
                  start = startglm,
                  iter,
                  epsilon,
                  report) %>% try(silent = F) -> fit2
  if("try-error" %in% class(fit2)){
    tauhat_glm  = NA; varhat_glm = NA
  }else{
    tauhat_glm = fit2$tau;
    varhat_glm = as.numeric(fit2$varhat)}
  re <- c(tauhat_gam_aic, tauhat_gam_bic, tauhat_glm,
          varhat_gam_aic, varhat_gam_bic, varhat_glm)
  names(re) <- c("tauhat_gam_aic", "tauhat_gam_bic", "tauhat_glm",
                 "varhat_gam_aic", "varhat_gam_bic", "varhat_glm")
  return(re)
}


simu_sense_new<- function(Ymodel,
                          Yname,
                          Xname,
                          Uname,
                          Tname,
                          emodel,
                          alpha0,
                          p01, p11,
                          s01, s11,
                          n,
                          v,
                          startgam,
                          startglm,
                          lambda,
                          gamma,
                          Kn,
                          p,
                          D,
                          ku ,
                          m,
                          iter,
                          epsilon,
                          report){
  dt1 <- Dtgeneration_new(Ymodel,
                          emodel,
                          alpha0,
                          p01, p11,
                          s01, s11 ,
                          n)  
  dt1 %>% filter(S==1) -> data1
  GAMmethod2_wrong(data1,
                   Yname,
                   Xname,
                   Uname,
                   Tname,
                   p01, p11 = 1, v, 1,
                   start = 0,
                   lambda,
                   gamma ,
                   Kn,
                   p,
                   D,
                   ku,
                   m,
                   iter,
                   epsilon,
                   report) %>% try(silent = F)-> fit1
  if("try-error" %in% class(fit1)){
    tauhat_gam_wr2 = NA
  }else{
    tauhat_gam_wr2 = fit1
  }
  GAMmethod2_wrong(data1,
                   Yname,
                   Xname,
                   Uname,
                   Tname,
                   p01, p11, v, s = 1,
                   start = 0,
                   lambda,
                   gamma ,
                   Kn,
                   p,
                   D,
                   ku,
                   m,
                   iter,
                   epsilon,
                   report) %>% try(silent = F)-> fit2
  if("try-error" %in% class(fit2)){
    tauhat_gam_wrs = NA
  }else{
    tauhat_gam_wrs = fit2
  }
  GAMmethod2_wrong(data1,
                   Yname,
                   Xname,
                   Uname,
                   Tname,
                   p01, p11 = 1, v, s = s11/s01,
                   start = startgam,
                   lambda,
                   gamma ,
                   Kn,
                   p,
                   D,
                   ku,
                   m,
                   iter,
                   epsilon,
                   report) %>% try(silent = F)-> fit3
  if("try-error" %in% class(fit3)){
    tauhat_gam_wrp = NA
  }else{
    tauhat_gam_wrp = fit3
  }
  GAMmethod2_wrong(data1,
                   Yname,
                   Xname,
                   Uname,
                   Tname,
                   p01, p11, v, s = s11/s01,
                   start = startgam,
                   lambda,
                   gamma ,
                   Kn,
                   p,
                   D,
                   ku,
                   m,
                   iter,
                   epsilon,
                   report) %>% try(silent = F)-> fit4
  if("try-error" %in% class(fit4)){
    tauhat_gam_ri = NA
  }else{
    tauhat_gam_ri = fit4
  }
  GLMmethod2_wrong(data1,
                   Yname,
                   Xname,
                   Uname,
                   Tname,
                   p01, p11=1, v, s=1,
                   start = 0,
                   iter,
                   epsilon,
                   report) %>% try(silent = F) -> fit5
  if("try-error" %in% class(fit5)){
    tauhat_glm_wr2  = NA; 
  }else{
    tauhat_glm_wr2 = fit5}
  GLMmethod2_wrong(data1,
                   Yname,
                   Xname,
                   Uname,
                   Tname,
                   p01, p11, v, s=1,
                   start = 0,
                   iter,
                   epsilon,
                   report) %>% try(silent = F) -> fit6
  if("try-error" %in% class(fit6)){
    tauhat_glm_wrs  = NA; 
  }else{
    tauhat_glm_wrs = fit6}  
  GLMmethod2_wrong(data1,
                   Yname,
                   Xname,
                   Uname,
                   Tname,
                   p01, p11=1, v, s=s11/s01,
                   start = startglm,
                   iter,
                   epsilon,
                   report) %>% try(silent = F) -> fit7
  if("try-error" %in% class(fit7)){
    tauhat_glm_wrp  = NA; 
  }else{
    tauhat_glm_wrp = fit7}  
  GLMmethod2_wrong(data1,
                   Yname,
                   Xname,
                   Uname,
                   Tname,
                   p01, p11, v, s=s11/s01,
                   start = startglm,
                   iter,
                   epsilon,
                   report) %>% try(silent = F) -> fit8
  if("try-error" %in% class(fit8)){
    tauhat_glm_ri  = NA; 
  }else{
    tauhat_glm_ri = fit8}    
  fit9 = summary(ate(Y_star ~ Tr | 1 | X1+X2+U, data=data1))[["asso"]][["coefmat"]]
  tau_iptw = -fit9[3,1]
  tau_van = IPTW_van(data1, v)
  tau_Di = tau_iptw/(p11-p01)
  return(c(tauhat_gam_wr2, tauhat_gam_wrs, tauhat_gam_wrp, tauhat_gam_ri,
           tauhat_glm_wr2, tauhat_glm_wrs, tauhat_glm_wrp, tauhat_glm_ri,
           tau_iptw , tau_Di, tau_van))
}
IPTW_van <- function(data1, v){
  data1$weight = ifelse(data1$Y_star, v, 1-v)
  glm(Tr~X1+X2+U, data = data1, family = binomial(), weights = weight) -> fit2
  predict(fit2, type = c("response")) -> data1$e_hat
  mean((data1$Tr/data1$e_hat - (1-data1$Tr)/(1-data1$e_hat))*
         data1$Y_star*data1$weight)/mean(data1$weight)  -> tau_van
  return(tau_van)
}

