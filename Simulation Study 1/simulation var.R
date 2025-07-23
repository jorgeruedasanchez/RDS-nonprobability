setwd(".../datos")

set.seed(1)
pob=read.table("fauxsycamore50.txt")
colnames(pob)=c("ips", "id","recruiter.id", "DEGREE", "y",
                "TONONDISEASE", "TODONDISEASE", "WAVE", "SEED", "WEIGHT")
N=nrow(pob)

total_y <- sum(pob$y)
prop_y <- total_y/nrow(pob)

est_RDS1=read.table("RDS-I estimates.txt")
est_RDS2=read.table("RDS-II estimates.txt")
est_RDSSS=read.table("RDS-SS estimates.txt")

vanal_ipw=vector()
vanal_mi=vector()
vanal_dr=vector()

#la elección del nr está en el main.

rep=1000
est=matrix(nrow=rep, ncol=6)
colnames(est) <- c("RDS-I", "RDS-II", "RDS-SS", "PSA1", "SM", "DR")

st_simul=Sys.time()
for(i in 1:rep){
  st_iter=Sys.time()
  set.seed(1+i)
  setwd(file.path(paste0("Sample ",i)))
  sv=read.table(file="muestra50repl.txt", header = FALSE, sep = "", dec = ".") #muestra non-prob (RDS)
  colnames(sv)=c("ips", "id","recruiter.id", "DEGREE", "y",
                 "TONONDISEASE", "TODONDISEASE", "WAVE", "SEED", "WEIGHT")
  
  setwd(".../datos")
  
  sr=pob[sample(nrow(pob), n_r, replace = F),]; #muestra prob (MAS)
  colnames(sr)=c("ips", "id","recruiter.id", "DEGREE", "y",
                 "TONONDISEASE", "TODONDISEASE", "WAVE", "SEED", "WEIGHT")
  
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
  
  
  est_PSA1 <- propest("y", sv, valliant_weights(propensidades$convenience))
  
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
  
  est[i,] <- c(est_RDS1[i,], est_RDS2[i,], est_RDSSS[i,], est_PSA1[[1]], est_SM, est_DR)
  
  # Var. analíticas
  ps_sv=propensidades$convenience
  ps_sr=propensidades$reference
  alg="glm"
  y_pred_sr=y_matching_sr
  y_pred_sv=y_matching_sv
  pesos="valliant"
  
  vanal_ipw[i]=var_anal_IPW(sv, sr, xr, xv, name_y, pi, ps_sv, ps_sr, pij_mat, pesos)
  vanal_mi[i]=var_anal_MI(sv, sr, xr, xv, name_y, pi, ps_sv, y_pred_sr, y_pred_sv, alg, pij_mat)
  an_dr=var_anal_DR(sv, sr, name_y, pi, ps_sv, ps_sr, y_pred_sv, y_pred_sr, pij_mat, pesos)
  vanal_dr[i]=an_dr$v1
  
  dfvar=data.frame(vanal_ipw, vanal_mi, vanal_dr)
  
  nombre_image = paste0("results_var.RData")
  save.image(file = nombre_image)

  print(paste0("iteración ",i, " de ", rep, " del algoritmo ", alg))
}

est=as.data.frame(est)

## Medidas de precisión

RBias=abs(apply(est,2,sum)/rep-prop_y)*(100/prop_y)
RSD=apply(est,2,sd)*(100/prop_y)
RMSE=sqrt(RBias^(2) + RSD^(2))

results=data.frame(RBias, RSD, RMSE)
