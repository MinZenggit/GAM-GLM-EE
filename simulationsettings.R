##########################
#### simulation setting!!
gets <- function(n, casesize, controlsize, p01, p11, v){
  # According the case-control size to calculate the selection probabilities for data generation.
  # v is the disease prevalence
  # p01 := P(Y*=1|Y=0); p11 := P(Y*=1|Y=1)
  # casesize/controlsize: the size of selected case and control samples from source population
  # n: the source population number.
  # s01 := P(S=1|Y*=0); s11 := P(S=1|Y*=1)
  v_star = (p11-p01)*v+p01
  s01 = controlsize/((1-v_star)*n)
  s11 = casesize/(v_star*n)
  return(c(s01, s11))
}

# disease prevalence
v.XS = 0.001
v.S = 0.005
v.L = 0.01
v.M = 0.05
v.H = 0.10
# Dtpara
emodel = 1
# GAMpara
Kn = 10
p = 3; D = 2; ku = 2; m = 2
lambda = c(5, 5); gamma = 0.01
iter = 200; epsilon = 0.000001
cores = 5; report = F
trials = 300
a1.XS = alpha0list[1]
a1.S = alpha0list[2]
a1.L = alpha0list[3]
a1.M = alpha0list[4]
a1.H = alpha0list[5]
a2.XS = alpha0list[6]
a2.S = alpha0list[7]
a2.L = alpha0list[8]
a2.M = alpha0list[9]
a2.H = alpha0list[10]
a3.XS = alpha0list[11]
a3.S = alpha0list[12]
a3.L = alpha0list[13]
a3.M = alpha0list[14]
a3.H = alpha0list[15]
a4.XS = alpha0list[16]
a4.S = alpha0list[17]
a4.L = alpha0list[18]
a4.M = alpha0list[19]
a4.H = alpha0list[20]
a5.XS = alpha0list[21]
a5.S = alpha0list[22]
a5.L = alpha0list[23]
a5.M = alpha0list[24]
a5.H = alpha0list[25]
a6.XS = alpha0list[26]
a6.S = alpha0list[27]
a6.L = alpha0list[28]
a6.M = alpha0list[29]
a6.H = alpha0list[30]
a7.XS = alpha0list[31]
a7.S = alpha0list[32]
a7.L = alpha0list[33]
a7.M = alpha0list[34]
a7.H = alpha0list[35]
a8.XS = alpha0list[36]
a8.S = alpha0list[37]
a8.L = alpha0list[38]
a8.M = alpha0list[39]
a8.H = alpha0list[40]
a9.XS = alpha0list[41]
a9.S = alpha0list[42]
a9.L = alpha0list[43]
a9.M = alpha0list[44]
a9.H = alpha0list[45]
setting =list()
for(i in 1:9){
  for(j in c("XS","S","L", "M", "H")){
    setting = rbind(setting, paste0("Ymodel=", i, ";v = v.", j, ";alpha0=a", i, ".", j ))
  }
}

print(setting)
startgam <- list()
for(i in 1:45){
  startgam[[i]] = c(rep(0, D*(Kn+p)), -2, -1, alpha0list[[i]])
}
startglm = list()
for(i in 1:45){
  startglm[[i]] = c(-1,0,-2, -1,alpha0list[i])
}   
