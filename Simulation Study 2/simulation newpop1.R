set.seed(1)

setwd("...")
pob=read.csv("population.csv")
N=nrow(pob)

total_y <- sum(pob$var1)
prop_y <- total_y/nrow(pob)

e_RDS2=read.table("RDS_II_estimates.txt")
e_RDSSS=read.table("RDS_SS_estimates.txt")

est_RDS2=vector()
est_RDSSS=vector()
est_IPW=vector()
est_SM=vector()
est_DR=vector()

#la elección del nr está en el main.

rep=1000
est=data.frame()

st_simul=Sys.time()
for(i in 1:rep){
  set.seed(1+i)
  sv=read.csv(file=paste0("rds_sample_", i, ".csv")) #muestra non-prob (RDS)
  sr=pob[sample(nrow(pob), n_r, replace = F),]; #muestra prob (MAS)

  pi=rep(n_r/N, n_r)
  sr$pi=pi
  name_y="y"
  pij_mat=Pkl.Hajek.s(pi)
  X="degree"
  xr=as.matrix(sr[, X])
  xv=as.matrix(sv[, X])
  
  source("script_sim.R")
  
  #PSA
  di <- 1/sr$pi
  propensidades <- propensities(sv, sr, "degree", algorithm = "glm",
                                model_weights = NULL,
                                smooth = FALSE, proc = NULL, trControl = trainControl(method = "none"))
  
  
  est_IPW[i] <- propest("var1", sv, valliant_weights(propensidades$convenience))[[1]]
  
  #SM
  y_matching_sr <- matching(sv %>% mutate(var1 = factor(var1, levels = c(0, 1), labels = c("X0", "X1"))),
                            sr %>% mutate(var1 = factor(var1, levels = c(0, 1), labels = c("X0", "X1"))), 
                            "degree", "var1", positive_label = "X1",algorithm = "glm")
  y_matching_sv <- matching(sv %>% mutate(var1 = factor(var1, levels = c(0, 1), labels = c("X0", "X1"))),
                            sv %>% mutate(var1 = factor(var1, levels = c(0, 1), labels = c("X0", "X1"))), 
                            "degree", "var1", positive_label = "X1",algorithm = "glm")
  
  est_SM[i]=sum(y_matching_sr*di)/nrow(pob)
  
  #DR
  est_DR[i] <- sum((sv$var1 - y_matching_sv)/propensidades$convenience)/nrow(pob) + est_SM[i]
  
  #RDS
  est_RDS2[i]=e_RDS2[i,]
  est_RDSSS[i]=e_RDSSS[i,]
  
  est <- data.frame(est_RDS2, est_RDSSS, est_IPW, est_SM, est_DR)

  
  save.image(file = paste0("results", escenario, ".RData"))
  print(paste0("iteración ",i, " de ", rep, " del escenario ", escenario))
}

colnames(est) <- c("RDS-II", "RDS-SS", "IPW", "SM", "DR")

## Medidas de precisión

RBias=abs(apply(est,2,sum)/rep-prop_y)*(100/prop_y)
RSD=apply(est,2,sd)*(100/prop_y)
RMSE=sqrt(RBias^(2) + RSD^(2))
results=data.frame(RBias, RSD, RMSE)