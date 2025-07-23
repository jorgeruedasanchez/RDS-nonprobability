setwd("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Escenario Inicial (estudio simul art)/datos")

set.seed(1)
pob=read.table("fauxsycamore50.txt")
colnames(pob)=c("ips", "id","recruiter.id", "DEGREE", "y",
                "TONONDISEASE", "TODONDISEASE", "WAVE", "SEED", "WEIGHT")

total_y <- sum(pob$y)
prop_y <- total_y/nrow(pob)


est_RDS1=read.table("RDS-I estimates.txt")
est_RDS2=read.table("RDS-II estimates.txt")
est_RDSSS=read.table("RDS-SS estimates.txt")

#la elección del nr está en el main.

rep=1000
est=matrix(nrow=rep, ncol=11)
colnames(est) <- c("Baseline", "RDS-I", "RDS-II", "RDS-SS", "PSA1", "PSA2", "PSA3", "PSA4", "SM", "DR", "KW")
p_psa1=vector()
p_psa2=vector()
p_psa3=vector()
p_psa4=vector()

for(i in 1:rep){
  set.seed(1+i)
  setwd(file.path(paste0("Sample ",i)))
  sv=read.table(file="muestra50repl.txt", header = FALSE, sep = "", dec = ".") #muestra non-prob (RDS)
  colnames(sv)=c("ips", "id","recruiter.id", "DEGREE", "y",
                 "TONONDISEASE", "TODONDISEASE", "WAVE", "SEED", "WEIGHT")
  
  setwd("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Escenario Inicial (estudio simul art)/datos")
  
  sr=pob[sample(nrow(pob), n_r, replace = F),]; #muestra prob (MAS)
  colnames(sr)=c("ips", "id","recruiter.id", "DEGREE", "y",
                 "TONONDISEASE", "TODONDISEASE", "WAVE", "SEED", "WEIGHT")
  sr$pi=nrow(sr)/nrow(pob)
  
  source("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Escenario Inicial (estudio simul art)/script_sim.R")
  
  #Baseline
  prop_y_sv <- sum(sv$y)/nrow(sv)
  
  #PSA
  di <- 1/sr$pi
  propensidades <- propensities(sv, sr, "DEGREE", algorithm = "glm",
                                model_weights = NULL,
                                smooth = FALSE, proc = NULL, trControl = trainControl(method = "none"))
  
  
  est_PSA1 <- propest("y", sv, valliant_weights(propensidades$convenience))
  est_PSA2 <- propest("y", sv, sc_weights(propensidades$convenience))
  est_PSA3 <- propest("y", sv, lee_weights(propensidades$convenience, propensidades$reference))
  est_PSA4 <- propest("y", sv, vd_weights(propensidades$convenience, propensidades$reference))
  
  #SM
  y_matching_sr <- matching(sv %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))),
                            sr %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))), 
                            "DEGREE", "y", positive_label = "X1",algorithm = "glm")
  y_matching_sv <- matching(sv %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))),
                            sv %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))), 
                            "DEGREE", "y", positive_label = "X1",algorithm = "glm")
  
  est_SM=sum(y_matching_sr*di)/nrow(pob); est_SM
  
  #DR
  est_DR <- sum((sv$y - y_matching_sv)/propensidades$convenience)/nrow(pob) + est_SM
  
  #KW
  wi_KW <- kw.wt(p_score.c = propensidades$convenience, p_score.s = propensidades$reference,
                 svy.wt = di)
  est_KW <- propest("y", sv, wi_KW$pswt)
  
  #Estimadores
  est[i,] <- c(prop_y_sv, est_RDS1[i,], est_RDS2[i,], est_RDSSS[i,], est_PSA1[[1]], est_PSA2[[1]], 
               est_PSA3[[1]], est_PSA4[[1]], est_SM, est_DR, est_KW[[1]])
  
  #Pesos
  p_psa1=c(p_psa1, valliant_weights(propensidades$convenience))
  p_psa2=c(p_psa2, sc_weights(propensidades$convenience))
  p_psa3=c(p_psa3, lee_weights(propensidades$convenience, propensidades$reference))
  p_psa4=c(p_psa4, vd_weights(propensidades$convenience, propensidades$reference))
  
  cat("Simulación nº", i, "con nr=", n_r, "\n")
}

est=as.data.frame(est)

## Medidas de precisión

RBias=abs(apply(est,2,sum)/rep-prop_y)*(100/prop_y)
RSD=apply(est,2,sd)*(100/prop_y)
RMSE=sqrt(RBias^(2) + RSD^(2))

results=data.frame(RBias, RSD, RMSE)

library(xlsx)
library(summarytools)
setwd("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Escenario Inicial (estudio simul art)/nr=1000 PSA")
pesos=cbind(p_psa1,p_psa2,p_psa3,p_psa4)
head(pesos)
## PSA2 es con make_preprocess_estimator y KW2 es con XGBoost.

descrip=descr(pesos); descrip
xlsx::write.xlsx(descrip, "pesos_summary.xlsx")

boxplot(pesos, main="Boxplot for weight comparison")
hist(pesos[,1], prob=T, main="Histogram for weight comparison", xlab = "IPW_weights")
lines(density(pesos[,1]),col="red",lwd=2)
lines(density(pesos[,2]),col="yellow",lwd=2)
lines(density(pesos[,3]),col="green",lwd=2)
lines(density(pesos[,4]),col="blue",lwd=2)

legend(x = "topright", legend = c("IPW/PSA1", "PSA2", "PSA3", "PSA4"), fill = c("red", "yellow", "green", "blue"))