##################
# This documents presented the codes for figure 3 in the main text of our paper.

library(Matrix)
library(extRC)
library(splines)
library(dplyr)
library(doParallel)
library(foreach)
library(rootSolve)
library(targeted)
rm(list = ls())
load("True_values_for_simulation/alpha0list.RData") # True alpha0 values
load("True_values_for_simulation/Ttau.RData") # True ate values
source("Dtgeneration_new.R")                  # data generation function
source("GAMAUTOGLM_new.R")                    # the main GAM,GLM methods function
source("GAMGLM_wrong.R")                      # the naive methods that ignore OMis or CC information.
source("Variance_estiamtion_new.R").          # variance estimation function for GAM/GLM methods 
source("Fisher_scoring_Algorithm2.R").        # optimization algorithm for GAM/GLM methods 
source("IC3.R")                               # Information criterion calculation function for GAM method
source("simu_new.R")                          # functions for simulation, including data generation and ATE estimation.
source("Resulttidy.R")                        # functions for presents the final simulation results
source("simulationsettings.R")                # simulation settings
set.seed(2023)
set = c(4, 5)

###########################
p01=0.00; p11 = 0.8;
casesize = controlsize =1000; n = 500000;
name = "1000-0.00-0.8test"
if(! dir.exists(paste0("simulation-results-figures/",name))){
  dir.create(paste0("simulation-results-figures/",name))
}
cores = 4; report = F
trials = 500
for(i in set){
  starttime <- Sys.time()
  eval(parse(text = setting[[i]]))
  gets(n, casesize, controlsize, p01 = 0, p11 = 1, v = v)-> a
  s01 = a[1]; s11 = a[2]
  cl<-makeCluster(cores)
  registerDoParallel(cl)
  result2 <- foreach(icount(trials),
                     .combine = rbind,
                     .errorhandling = "pass",
                     .packages = c("splines","Matrix","extRC","dplyr",
                                   "rootSolve", "targeted"))%dopar%simu_sense_new(
                                     Ymodel,
                                     Yname = "Y_star",
                                     Xname = c("X1", "X2"),
                                     Uname = "U",
                                     Tname = "Tr",
                                     emodel,
                                     alpha0,
                                     p01, p11,
                                     s01, s11,
                                     n,
                                     v,
                                     startgam[[i]],
                                     startglm[[i]],
                                     lambda,
                                     gamma,
                                     Kn,
                                     p,
                                     D,
                                     ku ,
                                     m,
                                     iter,
                                     epsilon,
                                     report = F)
  endtime <- Sys.time()
  stopCluster(cl)
  cat("setting ",i,": ")
  print(endtime-starttime)
  save(result2, file = paste0("simulation-results-figures/",name, "/", name, i,".rdata"))
}



###########################
p01=0.00; p11 = 0.6;
casesize = controlsize =1000; n = 500000;
name = "1000-0.00-0.6test"
if(! dir.exists(paste0("simulation-results-figures/",name))){
  dir.create(paste0("simulation-results-figures/",name))
}
cores = 4; report = F
trials = 500
for(i in set){
  starttime <- Sys.time()
  eval(parse(text = setting[[i]]))
  gets(n, casesize, controlsize, p01 = 0, p11 = 1, v = v)-> a
  s01 = a[1]; s11 = a[2]
  cl<-makeCluster(cores)
  registerDoParallel(cl)
  result2 <- foreach(icount(trials),
                     .combine = rbind,
                     .errorhandling = "pass",
                     .packages = c("splines","Matrix","extRC","dplyr",
                                   "rootSolve", "targeted"))%dopar%simu_sense_new(
                                     Ymodel,
                                     Yname = "Y_star",
                                     Xname = c("X1", "X2"),
                                     Uname = "U",
                                     Tname = "Tr",
                                     emodel,
                                     alpha0,
                                     p01, p11,
                                     s01, s11,
                                     n,
                                     v,
                                     startgam[[i]],
                                     startglm[[i]],
                                     lambda,
                                     gamma,
                                     Kn,
                                     p,
                                     D,
                                     ku ,
                                     m,
                                     iter,
                                     epsilon,
                                     report = F)
  endtime <- Sys.time()
  stopCluster(cl)
  cat("setting ",i,": ")
  print(endtime-starttime)
  save(result2, file = paste0("simulation-results-figures/",name, "/", name, i,".rdata"))
}







###############
## plot 
library(ggpubr)
library(tidyverse)
load("Ttau.RData")
plotresult_new <- function(result2, i, low, up, title){
  result2 <- result2 %>% subset(select = c(9, 1, 5, 10, 2, 6, 11, 3, 7, 4, 8))
  result2 %>% as.matrix() %>% as.numeric() -> a
  n = nrow(result2)
  class1 = rep(c("IPTW",   "GAM_wr2", "GLM_wr2", "IPTW_SY","GAM_wrs",
                 "GLM_wrs", "IPTW_Van","GAM_wrp", "GLM_wrp", "GAM_ri", "GLM_ri"), each = n)
  Method = c()
  for(j in 1:length(class1)){
    Method[j] = ifelse((str_detect(class1[j], "GAM")), "GAM-EE",
                       ifelse(str_detect(class1[j], "GLM"), "GLM-EE", 
                              ifelse(str_detect(class1[j], "SY"), "IPTW-SY", 
                                     ifelse(str_detect(class1[j], "Van"),"IPTW-Van", "IPTW"))))
  }
  Method =  factor(Method, levels = c("GAM-EE", "GLM-EE","IPTW","IPTW-SY","IPTW-Van"))
  class = c(rep(c( "None", "OMis","ODS"), each = 3*n), rep("ODS-OMis", 2*n))
  
  class = factor(class, levels = c("None", "OMis","ODS","ODS-OMis"))
  dat = data.frame(re = a, class = class, Method = Method)
  dat %>% ggplot( aes(x=class, y=re, fill=Method, alpha=Method)) + 
    geom_boxplot(varwidth = T, size = 0.2, outlier.size = 0.5) +
    scale_fill_manual(values=c("#69b3a2", "grey", "yellow", "red", "orange")) +
    scale_alpha_manual(values=c(0.8, 0.8, 0.5, 0.5, 0.5)) +theme_bw()+
    theme(legend.position = "top") + ylim(low, up) + ylab("ATE estimates")+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    ggtitle(title)+
    xlab("") + geom_hline(aes(yintercept=Ttau[i]), color = "red") -> picture
  return(picture)
}



name = "1000-0.00-0.8test"
i = 4
load(file = paste0("simulation-results-figures/", name,"/", name,i,".rdata"))
plotresult_new(result2, i, -0.6, 0, bquote(v == .(0.05) ~ "," ~ p[10] == .(0.2)))->p1
plot(p1)

name = "1000-0.00-0.8test"
i = 5
load(file = paste0("simulation-results-figures/", name,"/", name,i,".rdata"))
plotresult_new(result2, i, -0.6, 0, bquote(v == .(0.1) ~ "," ~ p[10] == .(0.2)))->p2
plot(p2)

name = "1000-0.00-0.6test"
i = 5
load(file = paste0("simulation-results-figures/", name,"/", name,i,".rdata"))
plotresult_new(result2, i, -0.67, 0, bquote(v == .(0.1) ~ "," ~ p[10] == .(0.4)))->p4
plot(p4)

name = "1000-0.00-0.6test"
i = 4
load(file = paste0("simulation-results-figures/", name,"/", name,i,".rdata"))
plotresult_new(result2, i, -0.67, 0, bquote(v == .(0.05) ~ "," ~ p[10] == .(0.4)))->p3
plot(p3)

# pdf(file = paste0("simulation-results-figures/", name, "/M1.pdf"), width = 20, height = 6, onefile = F)
# ggarrange(p2, p4,p1, p3,  nrow =1, ncol =4, common.legend  = T)
# dev.off()
