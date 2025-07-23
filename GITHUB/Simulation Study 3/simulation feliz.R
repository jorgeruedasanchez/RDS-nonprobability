#setwd("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Enc_indigenas/datos")
setwd("H:/Otros ordenadores/Mi portátil/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Enc_indigenas/datos") ## Desde portátil

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

total_y <- sum(pob_imput$feliz); med_y <- total_y/nrow(pob)
est_RDS1=read.table("RDS-I estimates y= feliz.txt")
est_RDS2=read.table("RDS-II estimates y= feliz.txt")
est_RDSSS=read.table("RDS-SS estimates y= feliz.txt")
rep=1000

est=matrix(nrow=rep, ncol=6)
colnames(est) <- c("RDS-I", "RDS-II", "RDS-SS", "PSA1", "SM", "DR")

#> mean(nv)
#[1] 161.862
#n_r=200
n_r=round(nrow(pob)*0.1); #nr=81

vanal_ipw=vector()
vanal_mi=vector()
vanal_dr=vector()

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

#setwd("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Enc_indigenas/datos")
setwd("H:/Otros ordenadores/Mi portátil/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Enc_indigenas/datos") ## Desde portátil

sr=pob_imput[sample(nrow(pob), n_r, replace = F),]; #muestra prob (MAS)
colnames(sr)=c("ips", "id", "recruiter.id", "degree", "sexo", "edad", "estadocivil", "cultura",
               "vestimenta", "instruccion", "internet", "idioma", "feliz", 
               "eleccion", "victima", "percepcion", "wave", "seed")

N=nrow(pob)
pi=rep(n_r/N, n_r)
sr$pi=pi
name_y="feliz"
pij_mat=Pkl.Hajek.s(pi)
X=c("degree", "sexo", "edad", "estadocivil", "cultura", "vestimenta", "instruccion",
    "internet", "idioma")
sr_x <- lapply(sr[,X], as.numeric)
sv_x <- lapply(sv[,X], as.numeric)

xr <- as.matrix(as.data.frame(sr_x))
xv <- as.matrix(as.data.frame(sv_x))


#source("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Enc_indigenas/vuelta al artículo/Revisión/script_sim.R")
source("H:/Otros ordenadores/Mi portátil/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Enc_indigenas/vuelta al artículo/Revisión/script_sim.R")

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

# Var. analíticas
ps_sv=propensidades$convenience
ps_sr=propensidades$reference
alg="glm"
y_pred_sr=y_matching_sr
y_pred_sv=y_matching_sv
pesos="valliant"

vanal_ipw[i]=var_anal_IPW(sv, sr, xr, xv, name_y, pi, ps_sv, ps_sr, pij_mat, pesos)
vanal_mi[i]=var_anal_MI(sv, sr, xr, xv, name_y, pi, ps_sv, y_pred_sr, y_pred_sv, alg, pij_mat)
an_dr=var_anal_DR(sv, sr, name_y, pi, ps_sv,ps_sr, y_pred_sv, y_pred_sr, pij_mat, pesos)
vanal_dr[i]=an_dr$v1

dfvar=data.frame(vanal_ipw, vanal_mi, vanal_dr)

nombre_image = paste0("resultsfeliz_rev.RData")
save.image(file = nombre_image)

outfile <- paste0("resultsfeliz_rev.RDataTmp")
if (file.exists(outfile)) {
  file.remove(outfile)
}
print(paste0("iteración ",i, " de ", rep, " del algoritmo ", alg))
}

est=as.data.frame(est)

## Medidas de precisión

RBias=abs(apply(est,2,sum)/rep-med_y)*(100/med_y)
RSD=apply(est,2,sd)*(100/med_y)
RMSE=sqrt(RBias^(2) + RSD^(2))

results_feliz=data.frame(RBias, RSD, RMSE)

## IC
nivel_confianza <- 0.95
alpha <- 1 - nivel_confianza
z <- qnorm(1 - alpha / 2)

CII_IPW=mean(est$PSA1-z*sqrt(dfvar$vanal_ipw))
CII_MI=mean(est$SM-z*sqrt(dfvar$vanal_mi))
CII_DR=mean(est$DR-z*sqrt(dfvar$vanal_dr))

CIS_IPW=mean(est$PSA1+z*sqrt(dfvar$vanal_ipw)) #IPW: 1.181569-1.301338
CIS_MI=mean(est$SM+z*sqrt(dfvar$vanal_mi)) #MI: 1.159242-1.318827
CIS_DR=mean(est$DR+z*sqrt(dfvar$vanal_dr)) #DR: 1.187808-1.290273
#Theta_mean=1.2432

