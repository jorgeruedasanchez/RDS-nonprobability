set.seed(1)
setwd("...")
pob=read.table("fauxsycamore50.txt")
colnames(pob)=c("ips", "id","recruiter.id", "DEGREE", "y",
                "TONONDISEASE", "TODONDISEASE", "WAVE", "SEED", "WEIGHT")
N=nrow(pob)

total_y <- sum(pob$y)
prop_y <- total_y/nrow(pob)

e_RDS2=read.table("RDS-II estimates.txt")
e_RDSSS=read.table("RDS-SS estimates.txt")

est_RDS2=vector()
est_RDSSS=vector()
est_IPW=vector()
est_SM=vector()
est_DR=vector()

rep=1000
est=data.frame()

st_simul=Sys.time()
for(i in 1:rep){
  set.seed(1+i)
  setwd(file.path(paste0("Sample ",i)))
  sv=read.table(file="muestra50repl.txt", header = FALSE, sep = "", dec = ".") #muestra non-prob (RDS)
  colnames(sv)=c("ips", "id","recruiter.id", "DEGREE", "y",
                 "TONONDISEASE", "TODONDISEASE", "WAVE", "SEED", "WEIGHT")
  setwd("...")
  
  sr=pob[sample(nrow(pob), n_r, replace = F),]; #muestra prob (MAS)

  pi=rep(n_r/N, n_r)
  sr$pi=pi
  name_y="y"
  pij_mat=Pkl.Hajek.s(pi)
  X="DEGREE"
  xr=as.matrix(sr[, X])
  xv=as.matrix(sv[, X])
  
  source("script_sim.R")

  #PSA
  di <- 1/sr$pi
  propensidades <- propensities(sv, sr, "DEGREE", algorithm = "glm",
                                model_weights = NULL,
                                smooth = FALSE, proc = NULL, trControl = trainControl(method = "none"))
  
  
  est_IPW[i] <- propest("y", sv, valliant_weights(propensidades$convenience))[[1]]
  
  #SM
  y_matching_sr <- matching(sv %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))),
                            sr %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))), 
                            "DEGREE", "y", positive_label = "X1",algorithm = "glm")
  y_matching_sv <- matching(sv %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))),
                            sv %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))), 
                            "DEGREE", "y", positive_label = "X1",algorithm = "glm")
  
  est_SM[i]=sum(y_matching_sr*di)/nrow(pob)
  
  #DR
  est_DR[i] <- sum((sv$y - y_matching_sv)/propensidades$convenience)/nrow(pob) + est_SM[i]
  
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