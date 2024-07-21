

getbias <- function(result, Ttau, d=3){
  a = apply(result,1,function(x)return(is.numeric(unlist(x))))
  convergrate <- mean(a)
  (result[a, ]) %>% unlist() %>% matrix(ncol = ncol(result))-> result 
  tauhat = mean(result[, 1])
  emsd <- sd(result[,1])
  sdhat = mean(sqrt(result[,2]))
  (result[,1]-Ttau)^2 %>% mean() %>% sqrt() -> RMSE
  Rbias = (tauhat-Ttau)/Ttau
  Lb <- result[,1]-sqrt(result[,2])*1.96
  Ub <- result[,1]+sqrt(result[,2])*1.96
  coverrate = mean(Ttau<Ub&Ttau>Lb)
  tdresult = c(tauhat,Rbias, coverrate, RMSE, convergrate)
  names(tdresult)=c("tauhat", "Rbias", "coverrate", "RMSE", "convergrate")
  return(round(tdresult, d))
}

getbias2 <- function(result, Ttau, d=3){
  convergrate <- 1-mean(is.na(result[,1]))
  na.omit(result) -> result
  tauhat = mean(result[, 1])
  emsd <- sd(result[,1])
  sdhat = mean(sqrt(result[,2]))
  (result[,1]-Ttau)^2 %>% mean() %>% sqrt() -> RMSE
  Rbias = (tauhat-Ttau)/Ttau
  Lb <- result[,1]-sqrt(result[,2])*1.96
  Ub <- result[,1]+sqrt(result[,2])*1.96
  coverrate = mean(Ttau<Ub&Ttau>Lb)
  tdresult = c(tauhat,Rbias, coverrate, RMSE, convergrate)
  names(tdresult)=c("tauhat", "Rbias", "coverrate", "RMSE", "convergrate")
  return(round(tdresult, d))
}



# # 
# # name = "11-3-1result"
# # name = "11-4-1result"
# # Ttau = c(Ttau1.L,Ttau1.M,Ttau1.H,
# #          Ttau2.L,Ttau2.M,Ttau2.H,
# #          Ttau3.L,Ttau3.M,Ttau3.H,
# #          Ttau4.L,Ttau4.M,Ttau4.H,
# #          Ttau5.L,Ttau5.M,Ttau5.H)
# # save(Ttau, file = "Ttau.RData")
# show <- c()
# for(i in 1:15){
#   load(file = paste0(name, "/settting", i,".RData"))
#   bias = c(getbias(result[,c(1, 2)], Ttau[i]), 
#            getbias(result[,c(3, 4)], Ttau[i]))
#   names(bias) = c("tauhat", "Rbias", "coverrate", "RMSE", "convergrate",
#                   "tauhat_glm", "Rbias_glm", "coverrate_glm", "RMSE_glm", "convergrate_glm")
#   show = rbind(show, bias)
# }
# rownames(show) <- paste0("Ymodel", rep(1:5, each = 3), c(".L", ".M", ".H"))
# 
# cbind(Ttau, show)  %>% round(digits = 3) -> show
# write.csv(show, file = paste0(name, ".csv"))
# 
# library(dplyr)
# load("~/Desktop/simu/11-7-4result5.RData")
# c(getbias(result[,c(1, 2)], Ttau7.M), 
#   getbias(result[,c(3, 4)], Ttau7.M))


