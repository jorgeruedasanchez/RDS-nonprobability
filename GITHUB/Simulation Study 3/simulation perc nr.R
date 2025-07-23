## Mismo script que SMML pero nr>nv
setwd("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Enc_indigenas/datos")
#setwd("H:/Otros ordenadores/Mi portátil/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Enc_indigenas/datos") ## Desde portátil

set.seed(1)
pob=read.table("indigenas.txt")
colnames(pob)=c("ips", "id", "recruiter.id", "degree", "sexo", "edad", "estadocivil", "cultura", "nacionalidad",
                "vestimenta", "instruccion", "propiedad", "internet", "idioma", "categoria", "valor", "feliz", 
                "eleccion", "victima", "percepcion", "wave", "seed")

# De las covariables nacionalidad(73), propiedad(6), internet(4), categoria(375), valor(375) tienen NA.
pob=subset(pob,select=-c(nacionalidad, categoria, valor, propiedad)) #Quitamos también propiedad ya que no sabemos las categorias de respuesta
pob$sexo=as.factor(pob$sexo)
pob$estadocivil=as.factor(pob$estadocivil)
pob$estadocivil[which(pob$estadocivil=="viudo/a")]="0"
pob$estadocivil=droplevels(pob$estadocivil)

pob$cultura=as.factor(pob$cultura)
pob$cultura[which(pob$cultura=="2")]="1"  # Agrupamos los niveles de respuesta para que no falle el matching
pob$cultura=droplevels(pob$cultura)

pob$vestimenta=as.factor(pob$vestimenta)
pob$instruccion=as.factor(pob$instruccion)
pob$internet=as.factor(pob$internet)
pob$idioma=as.factor(pob$idioma)

library(mice)
aa=mice(pob, seed=1); cat("Imputación poblacion \n")
pob_imput=complete(aa) #Imputar valores NA de la población

total_y <- sum(pob_imput$percepcion); med_y <- total_y/nrow(pob)
est_RDS2=read.table("RDS-II estimates y= feliz.txt")
est_RDSSS=read.table("RDS-SS estimates y= feliz.txt")

rep=1000

est=matrix(nrow=rep, ncol=10)
colnames(est) <- c("Baseline", "RDS-II", "RDS-SS", "PSA1", "PSA2", "PSA3", "PSA4", "SM", "DR", "KW")

#> mean(nv)
#[1] 161.862
n_r=200
#n_r=round(nrow(pob)*0.1); #nr=81

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
  sv$cultura[which(sv$cultura=="2")]="1"  # Agrupamos los niveles de respuesta para que no falle el matching
  sv$cultura=droplevels(sv$cultura)
  
  sv$vestimenta=as.factor(sv$vestimenta)
  sv$instruccion=as.factor(sv$instruccion)
  sv$internet=as.factor(sv$internet)
  sv$idioma=as.factor(sv$idioma)
  
  bb=mice(sv); cat("Imputación sv \n")
  sv=complete(bb, seed=1)
  
  setwd("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Enc_indigenas/datos")
  #setwd("H:/Otros ordenadores/Mi portátil/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Enc_indigenas/datos") ## Desde portátil
  
  sr=pob_imput[sample(nrow(pob), n_r, replace = F),]; #muestra prob (MAS)
  colnames(sr)=c("ips", "id", "recruiter.id", "degree", "sexo", "edad", "estadocivil", "cultura",
                 "vestimenta", "instruccion", "internet", "idioma", "feliz", 
                 "eleccion", "victima", "percepcion", "wave", "seed")
  
  sr$pi=nrow(sr)/nrow(pob)
  
  source("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Enc_indigenas/vuelta al artículo/script_sim_ML.R")
  #source("H:/Otros ordenadores/Mi portátil/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Enc_indigenas/vuelta al artículo/script_sim_ML.R")
  
  
  med_y_sv <- sum(sv$percepcion)/nrow(sv)
  
  ##################################################################################################################################
  ##  GLM ##
  ##########
  
  #PSA
  di <- 1/sr$pi
  propensidades <- propensities(sv, sr, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                          "internet", "idioma"),
                                algorithm = "glm",
                                model_weights = NULL,
                                smooth = FALSE, proc = NULL, trControl = trainControl(method = "none"))
  
  est_PSA1 <- meanest("percepcion", sv, valliant_weights(propensidades$convenience))
  est_PSA2 <- meanest("percepcion", sv, sc_weights(propensidades$convenience))
  est_PSA3 <- meanest("percepcion", sv, lee_weights(propensidades$convenience, propensidades$reference))
  est_PSA4 <- meanest("percepcion", sv, vd_weights(propensidades$convenience, propensidades$reference))
  
  #SM
  y_matching_sr <- matching_num(sv, sr, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                          "internet", "idioma"),
                                "percepcion", algorithm = "glm")
  y_matching_sv <- matching_num(sv, sv, c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
                                          "internet", "idioma"),
                                "percepcion", algorithm = "glm")
  
  est_SM=sum(y_matching_sr*di)/nrow(pob)
  
  #DR
  est_DR <- sum((sv$percepcion - y_matching_sv)/propensidades$convenience)/nrow(pob) + est_SM
  
  #KW
  wi_KW <- kw.wt(p_score.c = propensidades$convenience, p_score.s = propensidades$reference,
                 svy.wt = di)
  est_KW <- meanest("percepcion", sv, wi_KW$pswt)
  
  
  est[i,] <- c(med_y_sv, est_RDS2[i,], est_RDSSS[i,], est_PSA1[[1]], est_PSA2[[1]], est_PSA3[[1]], est_PSA4[[1]], est_SM, est_DR, est_KW[[1]])
  
  cat("Simulación nº", i, "\n")
}

est=as.data.frame(est)

## Medidas de precisión

RBias=abs(apply(est,2,sum)/rep-med_y)*(100/med_y)
RSD=apply(est,2,sd)*(100/med_y)
RMSE=sqrt(RBias^(2) + RSD^(2))

results_perc=data.frame(RBias, RSD, RMSE)