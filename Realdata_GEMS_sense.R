
###
# codes for real data analysis of GEMS dataset.
library(splines)
library(dplyr)
library(Matrix)
library(extRC)
library(doParallel)
library(foreach)
library(rootSolve)
library(targeted)
library(skimr)
rm(list = ls())
source("GAMAUTOGLM_new.R")
source("Variance_estiamtion_new.R")
source("Fisher_scoring_Algorithm2.R")
source("IC3.R")
source("Resulttidy.R")
load("ready.rdata")

dt_design = dt_ready
dt_design %>% select("Cryptos_ELISA", "Entamo_ELISA", "Giardia_ELISA", "Outcome",
                     "wealth_num", "wealth_cat", "washhand_before_nursing",
                     "washhand_after_defeca", "washhand_after_clean_kiddefe",
                     "Using_Soap_Yes", "Water_publictap", "Mali", "Gambia",
                     "Sex", "Age_group", "education", "weight_kg", "height_mean_cm",
                     "BMI_kgm2", "Age_month", "pet", "poultry", "cattle") -> dt_design





Vset = c(0.2, 0.22, 0.24, 0.26, 0.28, 0.3)
P01 = c(0.025, 0.020, 0.015, 0.010, 0.005)
P11 = sort(c(0.10, 0.15, 0.168, 0.2, 0.25, 0.3));
result1 <- list()
runningtime1 <- c()
result2 <- list()
runningtime2 <- c()

Kn = 5; p = 3; D = 2; ku = 8; m = 2;
lambdasetting = c(10, 10); gamma = 0.1;
iter = 1000; epsilon = 0.0000001;
k = 0


for(v in Vset){
  for(p11 in P11){
    for(p01 in P01){
      k = k+1
      # cat("case size =",case_size, "; control size =", control_size, "\n")
      startalpha_gam = -1.5
      startalpha_glm = 0
      cat("k=", k,"; p01 = ", p01, "; p11 = ", p11, "; v =", v,"\n")
      starttime = Sys.time()
      GAMmethod2_news(dt_design,
                      Yname = "Outcome",
                      Tname = "Cryptos_ELISA",
                      Uname = c("wealth_cat", "washhand_before_nursing", "Sex", "education", "Water_publictap", "Using_Soap_Yes","pet"),
                      Xname = c("BMI_kgm2", "Age_month"),
                      p01 = p01, p11 = p11, v = v,
                      start  = c(rep(0.0, D*(Kn+p)), 0, rep(0, ku-1),startalpha_gam),
                      lambda = lambdasetting,
                      gamma = gamma,
                      Kn = Kn,
                      p = p,
                      D = D,
                      ku = ku,
                      m = m,
                      iter = iter,
                      epsilon = epsilon,
                      report = F) %>% try(silent = T) -> fit_GAM
      result1[[k]] = fit_GAM
      endtime = Sys.time()
      runningtime1 = c(runningtime1, endtime - starttime)
      # cat("time1 = ", endtime - starttime, ";\ tauhat= ",fit_GAM$tauhat,"\n")
      starttime = Sys.time()
      GLMmethod2_news(dt_design,
                      Yname = "Outcome",
                      Tname = "Cryptos_ELISA",
                      Uname = c("wealth_cat", "washhand_before_nursing", "Sex", "education","Water_publictap","Using_Soap_Yes","pet"),
                      Xname = c("BMI_kgm2", "Age_month"),
                      p01 = p01, p11 = p11, v = v,
                      start = c(rep(0, ku+2), startalpha_glm),
                      iter = 200,
                      epsilon = 10^-7,
                      report = F) %>% try(silent = T) -> fit_GLM
      result2[[k]] = fit_GLM
      endtime = Sys.time()
      runningtime2 = c(runningtime2, endtime - starttime)
      # cat("time1 = ", endtime - starttime, ";\ tauhat= ",fit_GLM,"\n")
    }
  }
}



qn = matrix(0, nrow = length(Vset)*length(P11)*length(P01),
            ncol = 7)
qn %>%as.data.frame() -> qn
k = 0
for(v in Vset){
  for(p11 in P11){
    for(p01 in P01){
      k = k+1
      qn[k, 1] = v
      qn[k, 2] = p11
      qn[k, 3] = p01
      if(is.character(result1[[k]])==T){
        qn[k, 4] = NA
        qn[k, 5] = NA
      }else{
        qn[k, 4] = result1[[k]]$tauhat
        qn[k, 5] = result1[[k]]$varhat
      }
      if(is.character(result2[[k]])==T){
        qn[k, 6] = NA
        qn[k, 7] = NA
      }else{
        qn[k, 6] = result2[[k]]$tauhat
        qn[k, 7] = result2[[k]]$varhat
      }
      cat(k, ",")
    }}
}
colnames(qn) = c("v", "p11", "p01", "ate_a", "var_a", "ate_b", "var_b")
qn$low_a = qn$ate_a - 1.96*sqrt(qn$var_a)
qn$up_a = qn$ate_a + 1.96*sqrt(qn$var_a)
qn$low_b = qn$ate_b - 1.96*sqrt(qn$var_b)
qn$up_b = qn$ate_b + 1.96*sqrt(qn$var_b)

qn -> re
save(re, file = "GEMSresults/2-24-GEMS-2.rdata")


######################
##############################
library(ggpubr)
library(ggplot2)
library(dplyr)
library(latex2exp)
load("GEMSresults/2-24-GEMS-2.rdata")
################################
# re$ate_a <- ifelse(is.na(re$ate_a), re$ate_b, re$ate_a)
# re$low_a <- ifelse(is.na(re$low_a), re$low_b, re$low_a)
# re$up_a <- ifelse(is.na(re$up_a), re$up_b, re$up_a)

re$p10 = 1-re$p11
re %>% filter(p01!=0.012, p11!=0.168) -> re
re %>% filter(p01 == 0.005)  %>% select(p10, v, ate_a, low_a, up_a) -> dt
dt$p10 %>% as.character() -> dt$p10
ggplot(dt, aes(x = v, y = ate_a, group = p10, col = p10))+
  geom_line(position = position_dodge(0.005))+
  geom_point(position = position_dodge(0.005))+
  geom_errorbar(aes(ymin = low_a, ymax = up_a), width = 0.005,
                position = position_dodge(0.005)) +
  geom_hline(yintercept = 0, color = "black")+
  xlab("Disease prevelance")+
  ylab("ATE estimates") +ylim(c(0, 0.5))+
  scale_color_brewer(palette = "Blues")+
  theme_bw()+
  guides(col=guide_legend(title=bquote(p[.(sprintf("%02d", 10))])))+
  theme(legend.position = c(0.9, 0.25),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA))+
  ggtitle(bquote("GAM-EE:"~p[.(sprintf("%02d", 01))] == .(sprintf("%01.3f", 0.005)))) -> p1.0

re %>% filter(p01 == 0.01) %>% select(p10, v, ate_a, low_a, up_a) -> dt
dt$p10 %>% as.character() -> dt$p10
ggplot(dt, aes(x = v, y = ate_a, group = p10, col = p10))+
  geom_line(position = position_dodge(0.005))+
  geom_point(position = position_dodge(0.005))+
  geom_errorbar(aes(ymin = low_a, ymax = up_a), width = 0.005,
                position = position_dodge(0.005)) +
  geom_hline(yintercept = 0, color = "black")+
  xlab("Disease prevelance")+
  ylab("ATE estimates") +ylim(c(0, 0.5))+
  scale_color_brewer(palette = "Blues")+
  theme_bw()+
  guides(col=guide_legend(title=bquote(p[.(sprintf("%02d", 10))])))+
  theme(legend.position = c(0.9, 0.25),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA))+
  ggtitle(bquote("GAM-EE:"~p[.(sprintf("%02d", 01))] == .(sprintf("%01.3f", 0.01)))) -> p1.1

re %>% filter(p01 == 0.015) %>% select(p10, v, ate_a, low_a, up_a) -> dt
dt$p10 %>% as.character() -> dt$p10
ggplot(dt, aes(x = v, y = ate_a, group = p10, col = p10))+
  geom_line(position = position_dodge(0.005))+
  geom_point(position = position_dodge(0.005))+
  geom_errorbar(aes(ymin = low_a, ymax = up_a), width = 0.005,
                position = position_dodge(0.005)) +
  geom_hline(yintercept = 0, color = "black")+
  xlab("Disease prevelance")+
  ylab("ATE estimates") +ylim(c(0, 0.5))+
  scale_color_brewer(palette = "Blues")+
  theme_bw()+
  guides(col=guide_legend(title=bquote(p[.(sprintf("%02d", 10))])))+
  theme(legend.position = c(0.9, 0.25),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA))+
  ggtitle(bquote("GAM-EE:"~p[.(sprintf("%02d", 01))] == .(sprintf("%01.3f", 0.015)))) -> p1.2

re %>% filter(p01 == 0.020) %>% select(p10, v, ate_a, low_a, up_a) -> dt
dt$p10 %>% as.character() -> dt$p10
ggplot(dt, aes(x = v, y = ate_a, group = p10, col = p10))+
  geom_line(position = position_dodge(0.005))+
  geom_point(position = position_dodge(0.005))+
  geom_errorbar(aes(ymin = low_a, ymax = up_a), width = 0.005,
                position = position_dodge(0.005)) +
  geom_hline(yintercept = 0, color = "black")+
  xlab("Disease prevelance")+
  ylab("ATE estimates") + ylim(c(0, 0.5))+
  scale_color_brewer(palette = "Blues")+
  theme_bw()+
  guides(col=guide_legend(title=bquote(p[.(sprintf("%02d", 10))])))+
  theme(legend.position = c(0.9, 0.25),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA))+
  ggtitle(bquote("GAM-EE:"~p[.(sprintf("%02d", 01))] == .(sprintf("%01.3f", 0.020)))) -> p1.3


re %>% filter(p01 == 0.025) %>% select(p10, v, ate_a, low_a, up_a) -> dt
dt$p10 %>% as.character() -> dt$p10
ggplot(dt, aes(x = v, y = ate_a, group = p10, col = p10))+
  geom_line(position = position_dodge(0.005))+
  geom_point(position = position_dodge(0.005))+
  geom_errorbar(aes(ymin = low_a, ymax = up_a), width = 0.005,
                position = position_dodge(0.005)) +
  geom_hline(yintercept = 0, color = "black")+
  xlab("Disease prevelance")+
  ylab("ATE estimates") + ylim(c(0, 0.5))+
  scale_color_brewer(palette = "Blues")+
  theme_bw()+
  guides(col=guide_legend(title=bquote(p[.(sprintf("%02d", 10))])))+
  theme(legend.position = c(0.9, 0.25),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA))+
  ggtitle(bquote("GAM-EE:"~p[.(sprintf("%02d", 01))] == .(sprintf("%01.3f", 0.025)))) -> p1.4

pdf(file = "GEMSresults/GEMS-GAM-2-24-2.pdf", width = 12, height = 8)
ggarrange(p1.0, p1.1, p1.2, p1.3, p1.4, ncol = 3, nrow = 2,common.legend = T)
dev.off()


################################
re$p10 = 1-re$p11
re %>% filter(p01!=0.012, p11!=0.168) -> re
re %>% filter(p01 == 0.005)  %>% select(p10, v, ate_b, low_b, up_b) -> dt
dt$p10 %>% as.character() -> dt$p10
ggplot(dt, aes(x = v, y = ate_b, group = p10, col = p10))+
  geom_line(position = position_dodge(0.005))+
  geom_point(position = position_dodge(0.005))+
  geom_errorbar(aes(ymin = low_b, ymax = up_b), width = 0.005,
                position = position_dodge(0.005)) +
  geom_hline(yintercept = 0, color = "black")+
  xlab("Disease prevelance")+
  ylab("ATE estimates") +ylim(c(0, 0.6))+
  scale_color_brewer(palette = "Blues")+
  theme_bw()+
  guides(col=guide_legend(title=bquote(p[.(sprintf("%02d", 10))])))+
  theme(legend.position = c(0.9, 0.25),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA))+
  ggtitle(bquote("GLM-EE:"~p[.(sprintf("%02d", 01))] == .(sprintf("%01.3f", 0.005)))) -> p2.0

re %>% filter(p01 == 0.01) %>% select(p10, v, ate_b, low_b, up_b) -> dt
dt$p10 %>% as.character() -> dt$p10
ggplot(dt, aes(x = v, y = ate_b, group = p10, col = p10))+
  geom_line(position = position_dodge(0.005))+
  geom_point(position = position_dodge(0.005))+
  geom_errorbar(aes(ymin = low_b, ymax = up_b), width = 0.005,
                position = position_dodge(0.005)) +
  geom_hline(yintercept = 0, color = "black")+
  xlab("Disease prevelance")+
  ylab("ATE estimates") +ylim(c(0, 0.6))+
  scale_color_brewer(palette = "Blues")+
  theme_bw()+
  guides(col=guide_legend(title=bquote(p[.(sprintf("%02d", 10))])))+
  theme(legend.position = c(0.9, 0.25),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA))+
  ggtitle(bquote("GLM-EE:"~p[.(sprintf("%02d", 01))] == .(sprintf("%01.3f", 0.01)))) -> p2.1

re %>% filter(p01 == 0.015) %>% select(p10, v, ate_b, low_b, up_b) -> dt
dt$p10 %>% as.character() -> dt$p10
ggplot(dt, aes(x = v, y = ate_b, group = p10, col = p10))+
  geom_line(position = position_dodge(0.005))+
  geom_point(position = position_dodge(0.005))+
  geom_errorbar(aes(ymin = low_b, ymax = up_b), width = 0.005,
                position = position_dodge(0.005)) +
  geom_hline(yintercept = 0, color = "black")+
  xlab("Disease prevelance")+
  ylab("ATE estimates") +ylim(c(0, 0.6))+
  scale_color_brewer(palette = "Blues")+
  theme_bw()+
  guides(col=guide_legend(title=bquote(p[.(sprintf("%02d", 10))])))+
  theme(legend.position = c(0.9, 0.25),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA))+
  ggtitle(bquote("GLM-EE:"~p[.(sprintf("%02d", 01))] == .(sprintf("%01.3f", 0.015)))) -> p2.2

re %>% filter(p01 == 0.02) %>% select(p10, v, ate_b, low_b, up_b) -> dt
dt$p10 %>% as.character() -> dt$p10
ggplot(dt, aes(x = v, y = ate_b, group = p10, col = p10))+
  geom_line(position = position_dodge(0.005))+
  geom_point(position = position_dodge(0.005))+
  geom_errorbar(aes(ymin = low_b, ymax = up_b), width = 0.005,
                position = position_dodge(0.005)) +
  geom_hline(yintercept = 0, color = "black")+
  xlab("Disease prevelance")+
  ylab("ATE estimates") + ylim(c(0, 0.6))+
  scale_color_brewer(palette = "Blues")+
  theme_bw()+
  guides(col=guide_legend(title=bquote(p[.(sprintf("%02d", 10))])))+
  theme(legend.position = c(0.9, 0.25),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA))+
  ggtitle(bquote("GLM-EE:"~p[.(sprintf("%02d", 01))] == .(sprintf("%01.3f", 0.02)))) -> p2.3


re %>% filter(p01 == 0.025) %>% select(p10, v, ate_b, low_b, up_b) -> dt
dt$p10 %>% as.character() -> dt$p10
ggplot(dt, aes(x = v, y = ate_b, group = p10, col = p10))+
  geom_line(position = position_dodge(0.005))+
  geom_point(position = position_dodge(0.005))+
  geom_errorbar(aes(ymin = low_b, ymax = up_b), width = 0.005,
                position = position_dodge(0.005)) +
  geom_hline(yintercept = 0, color = "black")+
  xlab("Disease prevelance")+
  ylab("ATE estimates") + ylim(c(0, 0.6))+
  scale_color_brewer(palette = "Blues")+
  theme_bw()+
  guides(col=guide_legend(title=bquote(p[.(sprintf("%02d", 10))])))+
  theme(legend.position = c(0.9, 0.25),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA))+
  ggtitle(bquote("GLM-EE:"~p[.(sprintf("%02d", 01))] == .(sprintf("%01.3f", 0.025)))) -> p2.4


pdf(file = "GEMSresults/GEMS-GLM-2-24-2.pdf", width = 15, height = 4)
ggarrange(p2.0, p2.1, p2.2, p2.3, p2.4, ncol = 5, common.legend = T)
dev.off()

pdf(file = "GEMSresults/GEMS-GLM-GAM-2-25.pdf", width = 15, height = 8)
ggarrange(p1.0, p1.1, p1.2, p1.3, p1.4, 
          p2.0, p2.1, p2.2, p2.3, p2.4, ncol = 5, nrow = 2,common.legend = T)
dev.off()