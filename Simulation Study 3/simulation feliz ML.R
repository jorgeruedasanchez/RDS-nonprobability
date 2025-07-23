setwd(".../datos") 

set.seed(1)
pob=read.table("indigenas.txt")
colnames(pob)=c("ips", "id", "recruiter.id", "degree", "sexo", "edad", "estadocivil", "cultura", "nacionalidad",
                "vestimenta", "instruccion", "propiedad", "internet", "idioma", "categoria", "valor", "feliz", 
                "eleccion", "victima", "percepcion", "wave", "seed")

pob=subset(pob,select=-c(nacionalidad, categoria, valor, propiedad))
pob$sexo=as.factor(pob$sexo)
pob$estadocivil=as.factor(pob$estadocivil)
pob$estadocivil[which(pob$estadocivil=="viudo/a")]="0"
pob$estadocivil=droplevels(pob$estadocivil)

pob$cultura=as.factor(pob$cultura)
pob$cultura[which(pob$cultura=="2")]="1"  
pob$cultura=droplevels(pob$cultura)

pob$vestimenta=as.factor(pob$vestimenta)
pob$instruccion=as.factor(pob$instruccion)
pob$internet=as.factor(pob$internet)
pob$idioma=as.factor(pob$idioma)

library(mice)
aa=mice(pob, seed=1)
pob_imput=complete(aa) 

total_y <- sum(pob_imput$feliz); med_y <- total_y/nrow(pob)
est_RDS1=read.table("RDS-I estimates y= feliz.txt")
est_RDS2=read.table("RDS-II estimates y= feliz.txt")
est_RDSSS=read.table("RDS-SS estimates y= feliz.txt")
rep=1000

est=matrix(nrow=rep, ncol=6)
colnames(est) <- c("RDS-I", "RDS-II", "RDS-SS", "PSA1", "SM", "DR")

est_knn=matrix(nrow=rep, ncol=3)
colnames(est_knn) <- c("PSA1", "SM", "DR")

est_gbm=matrix(nrow=rep, ncol=3)
colnames(est_gbm) <- c("PSA1", "SM", "DR")

est_crf=matrix(nrow=rep, ncol=3)
colnames(est_crf) <- c("PSA1", "SM", "DR")


n_r=round(nrow(pob)*0.1); #nr=81

for(i in 1:rep){
set.seed(1+i)
setwd(file.path(paste0("Sample ",i)))
sv=read.table(file="muestra50repl.txt", header = FALSE, sep = "", dec = ".") #muestra non-prob (RDS)
colnames(sv)=c("ips", "id", "recruiter.id", "degree", "sexo", "edad", "estadocivil", "cultura", "nacionalidad",
               "vestimenta", "instruccion", "propiedad", "internet", "idioma", "categoria", "valor", "feliz", 
               "eleccion", "victima", "percepcion", "wave", "seed")

sv=subset(sv,select=-c(nacionalidad, categoria, valor, propiedad))
sv$sexo=as.factor(sv$sexo)
sv$estadocivil=as.factor(sv$estadocivil)
sv$estadocivil[which(sv$estadocivil=="viudo/a")]=0
sv$estadocivil=droplevels(sv$estadocivil)

sv$cultura=as.factor(sv$cultura)
sv$cultura[which(sv$cultura=="2")]="1"  
sv$cultura=droplevels(sv$cultura)

sv$vestimenta=as.factor(sv$vestimenta)
sv$instruccion=as.factor(sv$instruccion)
sv$internet=as.factor(sv$internet)
sv$idioma=as.factor(sv$idioma)

bb=mice(sv); cat("Imputación sv \n")
sv=complete(bb, seed=1)


setwd(".../datos") 

sr=pob_imput[sample(nrow(pob), n_r, replace = F),]; #muestra prob (MAS)
colnames(sr)=c("ips", "id", "recruiter.id", "degree", "sexo", "edad", "estadocivil", "cultura",
               "vestimenta", "instruccion", "internet", "idioma", "feliz", 
               "eleccion", "victima", "percepcion", "wave", "seed")

sr$pi=nrow(sr)/nrow(pob)

source("script_sim_ML.R")

##################################################################################################################################
##  GLM ##
##########
med_y_sv <- sum(sv$feliz)/nrow(sv)

#PSA
di <- 1/sr$pi
propensidades <- propensities(sv, sr, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                        "internet", "idioma"),
                              algorithm = "glm",
                              model_weights = NULL,
                              smooth = FALSE, proc = NULL, trControl = trainControl(method = "none"))

est_PSA1 <- meanest("feliz", sv, valliant_weights(propensidades$convenience))

#SM
y_matching_sr <- matching_num(sv, sr, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                        "internet", "idioma"),
                              "feliz", algorithm = "glm")
y_matching_sv <- matching_num(sv, sv, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                        "internet", "idioma"),
                              "feliz", algorithm = "glm")

est_SM=sum(y_matching_sr*di)/nrow(pob)

#DR
est_DR <- sum((sv$feliz - y_matching_sv)/propensidades$convenience)/nrow(pob) + est_SM


est[i,] <- c(est_RDS1[i,], est_RDS2[i,], est_RDSSS[i,], est_PSA1[[1]], est_SM, est_DR)

##################################################################################################################################
## K-vecinos ##
###############

#PSA
di <- 1/sr$pi
propensidades_knn <- propensities(sv, sr, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                            "internet", "idioma"),
                                  algorithm = "knn",
                                  model_weights = NULL,
                                  smooth = FALSE, proc = NULL, trControl = trainControl(method = "none"))

est_PSA1_knn <- meanest("feliz", sv, valliant_weights(propensidades_knn$convenience))

#SM
y_matching_sr_knn <- matching_num(sv, sr, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                            "internet", "idioma"),
                                  "feliz", algorithm = "knn")
y_matching_sv_knn <- matching_num(sv, sv, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                            "internet", "idioma"),
                                  "feliz", algorithm = "knn")

est_SM_knn=sum(y_matching_sr_knn*di)/nrow(pob)

#DR
est_DR_knn <- sum((sv$feliz - y_matching_sv_knn)/propensidades_knn$convenience)/nrow(pob) + est_SM_knn


est_knn[i,] <- c(est_PSA1_knn[[1]], est_SM_knn, est_DR_knn)


##################################################################################################################################
## GBM ##
#########

#PSA
di <- 1/sr$pi
propensidades_gbm <- propensities(sv, sr, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                            "internet", "idioma"),
                                  algorithm = "gbm",
                                  model_weights = NULL,
                                  smooth = FALSE, proc = NULL, trControl = trainControl(method = "none"))

est_PSA1_gbm <- meanest("feliz", sv, valliant_weights(propensidades_gbm$convenience))

#SM
y_matching_sr_gbm <- matching_num(sv, sr, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                            "internet", "idioma"),
                                  "feliz", algorithm = "gbm")
y_matching_sv_gbm <- matching_num(sv, sv, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                            "internet", "idioma"),
                                  "feliz", algorithm = "gbm")

est_SM_gbm=sum(y_matching_sr_gbm*di)/nrow(pob)

#DR
est_DR_gbm <- sum((sv$feliz - y_matching_sv_gbm)/propensidades_gbm$convenience)/nrow(pob) + est_SM_gbm


est_gbm[i,] <- c(est_PSA1_gbm[[1]], est_SM_gbm, est_DR_gbm)

##################################################################################################################################
## Cond. Random forest ##
#########################

#PSA
di <- 1/sr$pi
propensidades_crf <- propensities(sv, sr, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                             "internet", "idioma"),
                                   algorithm = "cforest",
                                   model_weights = NULL,
                                   smooth = FALSE, proc = NULL, trControl = trainControl(method = "none"))

est_PSA1_crf <- meanest("feliz", sv, valliant_weights(propensidades_crf$convenience))

#SM
y_matching_sr_crf <- matching_num(sv, sr, c("sexo", "estadocivil", "cultura", "vestimenta", "instruccion",
                                             "internet", "idioma"),
                                   "feliz", algorithm = "cforest")
y_matching_sv_crf <- matching_num(sv, sv, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                             "internet", "idioma"),
                                   "feliz", algorithm = "cforest")

est_SM_crf=sum(y_matching_sr_crf*di)/nrow(pob)

#DR
est_DR_crf <- sum((sv$feliz - y_matching_sv_crf)/propensidades_crf$convenience)/nrow(pob) + est_SM_crf


est_crf[i,] <- c(est_PSA1_crf[[1]], est_SM_crf, est_DR_crf)


cat("Simulación nº", i, "\n")
}

est=as.data.frame(est)
est_knn=as.data.frame(est_knn)
est_gbm=as.data.frame(est_gbm)
est_crf=as.data.frame(est_crf)


## Medidas de precisión

RBias=abs(apply(est,2,sum)/rep-med_y)*(100/med_y)
RSD=apply(est,2,sd)*(100/med_y)
RMSE=sqrt(RBias^(2) + RSD^(2))

results_feliz=data.frame(RBias, RSD, RMSE)

#knn
RBias_knn=abs(apply(est_knn,2,sum)/rep-med_y)*(100/med_y)
RSD_knn=apply(est_knn,2,sd)*(100/med_y)
RMSE_knn=sqrt(RBias_knn^(2) + RSD_knn^(2))

knn_feliz=data.frame(RBias_knn, RSD_knn, RMSE_knn)

#gbm
RBias_gbm=abs(apply(est_gbm,2,sum)/rep-med_y)*(100/med_y)
RSD_gbm=apply(est_gbm,2,sd)*(100/med_y)
RMSE_gbm=sqrt(RBias_gbm^(2) + RSD_gbm^(2))

gbm_feliz=data.frame(RBias_gbm, RSD_gbm, RMSE_gbm)

#CRF
RBias_crf=abs(apply(est_crf,2,sum)/rep-med_y)*(100/med_y)
RSD_crf=apply(est_crf,2,sd)*(100/med_y)
RMSE_crf=sqrt(RBias_crf^(2) + RSD_crf^(2))

crf_feliz=data.frame(RBias_crf, RSD_crf, RMSE_crf)