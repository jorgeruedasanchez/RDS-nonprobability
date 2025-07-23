library(tidyverse)
library(lme4)
library(caret)
library(NonProbEst)
library(caret)
library(sampling)
library(survey)
library(FSelector)
library(KWML)
library(xlsx)

#### SAMPLING ####

muestreo.lahiri <- function(relsize, n){
  C <- sum(tail(sort(relsize), n))
  cumple <- 0
  while(cumple == 0){
    num_e <- runif(1, 0, C)
    muestra_conr <- sample(length(relsize), n, replace = F)
    if(sum(relsize[muestra_conr]) >= num_e) cumple <- 1
  }
  return(muestra_conr)
}

muestreo_cohorte <- function(n, m){
  muestra_fase1 <- muestreo.lahiri(ps_fase1_cohorte, n)#which(sampling::UPmidzuno(ps_fase1*n) >= 1)
  muestra_fase2 <- list()
  for(i in 1:length(muestra_fase1)){
    muestra_temp <- muestreo.lahiri(ps_fase2_cohorte[[muestra_fase1[i]]], m) #which(sampling::UPmidzuno(ps_fase2[[muestra_fase1[i]]]*m) >= 1)
    muestra_fase2[[i]] <- clusters[[muestra_fase1[i]]][muestra_temp,]
  }
  muestra_final <- do.call(rbind.data.frame, muestra_fase2)
  return(muestra_final)
}

muestreo_prob <- function(n, m){
  muestra_fase1 <- muestreo.lahiri(ps_fase1_prob, n)#which(sampling::UPmidzuno(ps_fase1*n) >= 1)
  muestra_fase2 <- list()
  for(i in 1:length(muestra_fase1)){
    muestra_temp <- muestreo.lahiri(ps_fase2_prob[[muestra_fase1[i]]], m) #which(sampling::UPmidzuno(ps_fase2[[muestra_fase1[i]]]*m) >= 1)
    muestra_fase2[[i]] <- clusters[[muestra_fase1[i]]][muestra_temp,]
  }
  muestra_final <- do.call(rbind.data.frame, muestra_fase2)
  return(muestra_final)
}

propest <- function(variable, s_v, wi_i) svyciprop(as.formula(paste0("~I(",variable,"==1)")), 
                                                   design = trimWeights(svydesign(~0, data = s_v, weights = wi_i), lower = 0), 
                                                   method = "mean")

matching <- function (convenience_sample, reference_sample, covariates, estimated_var, 
                      positive_label = NULL, algorithm = "glm", proc = NULL,
                      trControl = trainControl(classProbs = TRUE, method = "none"), ...) 
{
  data = convenience_sample[, covariates, drop = FALSE]
  values = convenience_sample[, estimated_var]
  test = reference_sample[, covariates, drop = FALSE]
  model = train(data, values, algorithm, preProcess = proc, 
                trControl = trControl, ...)
  if (is.null(positive_label)) 
    return(predict(model, test))
  else return(predict(model, test, type = "prob")[, positive_label])
}

model_todos <- function(sample_data, weights, full_data, covariates, estimated_var, 
                        estimate_mean = FALSE, positive_label = NULL, algorithm = "glm", 
                        proc = NULL, ...) 
{
  if (length(weights) == 1) 
    weights = rep(weights, nrow(sample_data))
  known_values = sample_data[, estimated_var]
  if (!is.null(positive_label)) 
    known_values = known_values == positive_label
  all_data = rbind(sample_data[, covariates, drop = FALSE], 
                   full_data[, covariates, drop = FALSE])
  all_predicted_values = matching(sample_data, all_data, covariates, 
                              estimated_var, positive_label, algorithm = algorithm, 
                              proc = proc, ...)
  known_predicted_values = all_predicted_values[1:length(known_values)]
  
  # Model-based
  predicted_values = all_predicted_values[(length(known_values) + 
                                             1):length(all_predicted_values)]
  total_mb = sum(known_values, predicted_values, -known_predicted_values)
  
  # Model-assisted
  total_ma = sum(all_predicted_values, -known_predicted_values) + 
    sum((known_values - known_predicted_values) * weights)
  
  # Model-calibrated
  predicted_total = sum(all_predicted_values, -known_predicted_values)
  final_weights = calib(known_predicted_values, weights, predicted_total, 
                        method = "linear") * weights
  total_mc = sum(known_values * final_weights)
  
  if (estimate_mean) 
    return(c(total_mb, total_ma, total_mc)/nrow(full_data))
  else return(c(total_mb, total_ma, total_mc))
}

model_hardsoft <- function(sample_data, weights, full_data, X1, X2, estimated_var, 
                           estimate_mean = FALSE, positive_label = NULL, 
                           proc = NULL, ...) 
{
  if (length(weights) == 1) 
    weights = rep(weights, nrow(sample_data))
  covariates <- c(X1, X2)
  familia <- "gaussian"
  known_values = sample_data[, estimated_var]
  if (!is.null(positive_label)){
    known_values = known_values == positive_label
    familia <- "binomial"
  }
  all_data = rbind(sample_data[, covariates, drop = FALSE], 
                   full_data[, covariates, drop = FALSE])
  formula_mlr <- paste0(estimated_var, " ~ ", paste(X1, collapse = " + "), 
                        " + (1|", paste(X2, collapse = ") + (1|"), ")")
  all_predicted_values = predict(glmer(formula_mlr, data = sample_data, family = familia), all_data)
  known_predicted_values = all_predicted_values[1:length(known_values)]
  
  # Model-based
  predicted_values = all_predicted_values[(length(known_values) + 
                                             1):length(all_predicted_values)]
  total_mb = sum(known_values, predicted_values, -known_predicted_values)
  
  # Model-assisted
  total_ma = sum(all_predicted_values, -known_predicted_values) + 
    sum((known_values - known_predicted_values) * weights)
  
  # Model-calibrated
  predicted_total = sum(all_predicted_values, -known_predicted_values)
  final_weights = calib(known_predicted_values, weights, predicted_total, 
                        method = "linear") * weights
  total_mc = sum(known_values * final_weights)
  
  if (estimate_mean) 
    return(c(total_mb, total_ma, total_mc)/nrow(full_data))
  else return(c(total_mb, total_ma, total_mc))
}

propensities <- function (convenience_sample, reference_sample, covariates, algorithm = "glm", model_weights = NULL,
                          smooth = FALSE, proc = NULL, trControl = trainControl(classProbs = TRUE), 
                          ...) 
{
  n_convenience = nrow(convenience_sample)
  n_reference = nrow(reference_sample)
  data = rbind(convenience_sample[, covariates, drop = FALSE], 
               reference_sample[, covariates, drop = FALSE])
  labels = append(rep(1, n_convenience), rep(0, n_reference))
  if(is.null(model_weights))  model_weights = append(rep(1, n_convenience), 
                                                     rep(n_convenience/n_reference, n_reference))
  trControl$classProbs = TRUE
  model = train(data, factor(labels, levels = c(1, 0), labels = c("Positive", 
                                                                  "Negative")), algorithm, weights = model_weights, preProcess = proc, 
                trControl = trControl, ...)
  probabilities = predict(model, data, type = "prob")$Positive
  if (smooth) 
    probabilities = (1000 * probabilities + 0.5)/1001
  list(convenience = probabilities[1:n_convenience], reference = probabilities[(n_convenience + 
                                                                                  1):length(probabilities)])
}

papp <- function (convenience_sample, reference_sample, covariates, weights, 
                  algorithm = "glm", proc = NULL, ...) 
{
  test = convenience_sample[, covariates, drop = FALSE]
  values = weights
  data = reference_sample[, covariates, drop = FALSE]
  model = train(data, values, algorithm, preProcess = proc, 
                trControl = trainControl(classProbs = TRUE, method = "none"), ...)
  return(predict(model, test))
}

psaplusmatch <- function(muestra_cohorte, muestra_prob, pesos, di, X, algoritmo = "glm",
                             trControl = trainControl(classProbs = TRUE, method = "none")){
  npesos <- nrow(muestra_cohorte)*pesos/sum(pesos)
  y_matching_psa <- matching(muestra_cohorte %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))),
                             muestra_prob %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))), 
                             X, "y", positive_label = "X1", algorithm = algoritmo, weights = npesos,
                             trControl = trControl)
  return(c(sum(y_matching_psa * di)/3000000, sum(y_matching_psa * di)/sum(di)))
}

#### VARIANZAS ANALÍTICAS ####

sumatorio_ipw=function(j, i, pij_mat, ps_sr, xr) {
  sum_D=((pij_mat[i,j]-(pij_mat[i,i]*pij_mat[j,j]))/pij_mat[i,j]) * (ps_sr[i]*ps_sr[j]/(pij_mat[i,i]*pij_mat[j,j])) 
  return(sum_D*((xr[i,])%*%t(xr[j,]))) #esto devuelve una matriz ncol*ncol 5*5
}

var_anal_IPW <- function (sv, sr, xr, xv, name_y, pi, ps_sv, ps_sr, pij_mat=NULL, pesos){
  clusterExport(cluster, ls(), envir=environment())
  #sv y sr deben ser data.frames
  #xr y xv deben ser matrices de variables indicadoras (columnas individuos, filas variables)  
  
  nr=nrow(sr)
  nv=nrow(sv)
  y=sv[, name_y]
  di=1/pi
  Nr=sum(di)
  
  if(is.null(pij_mat)==TRUE){
    pij_mat=Pkl.Hajek.s(pi)
  }
  if(pesos=="valliant")
    wi_ipw=valliant_weights(ps_sv) #CAMBIAR
  if(pesos=="sc")
    wi_ipw=sc_weights(ps_sv) #CAMBIAR
  if(pesos=="kernel")
    wi_ipw=kw.wt(p_score.c =ps_sv, p_score.s =ps_sr, svy.wt = N/n_A)$pswt #CAMBIAR 
  
  est_ipw= sum(y*wi_ipw)/sum(wi_ipw)
  Nv=sum(wi_ipw)
  
  ## D
  D=0
  sum_D=list()
  for(i in 1:nr){
    sum_D[[i]] <-parSapply(cluster, 1:nr,  sumatorio_ipw, i, pij_mat, ps_sr, xr) 
  }
  if (is.null(sum_D[[1]]) || !is.matrix(sum_D[[1]])) {
    # Si los elementos de sum_D son vectores (ej. de sumatorio_ipw devolviendo un solo número)
    message("Los elementos de sum_D son vectores; usando sum() para la agregación.")
    suma_filas <- lapply(sum_D, sum)
  } else {
    # Si los elementos de sum_D son matrices (ej. de sumatorio_ipw devolviendo un vector de longitud > 1)
    message("Los elementos de sum_D son matrices; usando rowSums() para la agregación.")
    suma_filas <- lapply(sum_D, rowSums)
  }
  
  #suma_filas <- lapply(sum_D, rowSums)#sumo por filas los elementos de la lista (matrices) obteniendo una lista de 15 matrices de 25*1 
  suma_por_posicion <- Reduce(`+`, suma_filas)#de las 15 matrices voy sumando los elementos por posicion
  
  D <- matrix(suma_por_posicion, nrow = ncol(xr), ncol = ncol(xr), byrow = TRUE) #lo reestructuro como una matriz 5*5 (columnas de xr var aux)
  
  #D_def=D/N^2
  D_def=D/Nr^2
  
  ## b2
  num_b2=0
  den_b2=0
  for(i in 1:nv){
    num_b2=num_b2 + (1/ps_sv[i]-1)*(y[i]-est_ipw)*t(xv[i,])
  }  
  for(i in 1:nr){
    den_b2=den_b2 + di[i]*(1-ps_sr[i])*ps_sr[i]*(xr[i,] %*% t(xr[i,]))
  }
  b2_t=num_b2 %*% solve(den_b2) 
  
  ## Fórmula var_ipw
  v_ipw_2=b2_t%*%D_def%*%t(b2_t); v_ipw_2
  
  v_ipw_1=0
  for(i in 1:nv){
    v_ipw_1=v_ipw_1 + (((y[i]-est_ipw)/ps_sv[i])-b2_t%*%xv[i,])^2 * ((1-ps_sv[i]))
  }
  
  #v_ipw_ED=(v_ipw_1/N^2)+v_ipw_2
  v_ipw=(v_ipw_1/Nv^2)+v_ipw_2
  
  return(v_ipw) #0.01284937
}

sumatorio_mi_dr=function(j,i,pij_mat, di, y_pred_sr) {
  suma=((pij_mat[i,j]-(pij_mat[i,i]*pij_mat[j,j]))/pij_mat[i,j]) * (di[i]*y_pred_sr[i]*di[j]*y_pred_sr[j])
  return(suma)
  
}

var_anal_MI <- function (sv, sr, xr, xv, name_y, pi, ps_sv, y_pred_sr,  y_pred_sv=NULL, alg, pij_mat=NULL) {
  
  clusterExport(cluster, ls(), envir=environment()) 
  #sv y sr deben ser data.frames
  #xr y xv deben ser matrices de variables indicadoras (columnas individuos, filas variables)  
  
  nr=nrow(sr)
  nv=nrow(sv)
  y=sv[, name_y]
  di=1/pi
  N=sum(di)
  
  
  if(is.null(pij_mat)==TRUE){
    pij_mat=Pkl.Hajek.s(pi)
  }
  
  if(is.null(y_pred_sv)==TRUE){
    if(length(unique(sv[,`name_y`]))==2){
      y_pred_sv <- matching(sv %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))),
                            sv, X, name_y, positive_label = "X1", algorithm = alg)
    } else {
      y_pred_sv <- matching(sv, sv, X, name_y, algorithm = alg)
    }
  }
  
  ## c
  den_c=0
  for(i in 1:nv){
    den_c=den_c+ t(xv[i,])%*%xv[i,]
  }
  den_c=(1/nv)*den_c
  
  num_c=0
  for(i in 1:nr){
    num_c=num_c + di[i]*xr[i,]
  }
  c=(1/N)*(num_c%*%(den_c^(-1)))
  
  
  ## e (residuos)
  e=y-y_pred_sv
  
  ## V1
  v1=0
  for(i in 1:nv){
    v1=v1 + (e[i]^2)*(t(xv[i,])%*%c)^2
  }
  v_mi_1=v1/(nv^2)
  
  
  ## V2
  v2=0
  sum_v2=matrix(NA, nr, nr)
  for(i in 1:nr){
    sum_v2[i,]=parSapply(cluster, 1:nr,   sumatorio_mi_dr, i, pij_mat, di, y_pred_sr) #matriz nr*nr
  }
  v2=sum(sum_v2) 
  v_mi_2=v2/N^2
  
  ## vAR FINAL
  v_mi=v_mi_2+v_mi_1
  
  return(v_mi)  #0.007171927
}


var_anal_DR <- function (sv, sr, name_y, pi, ps_sv, ps_sr, y_pred_sv, y_pred_sr, pij_mat=NULL, pesos) {
  clusterExport(cluster, ls(), envir=environment()) 
  #sv y sr deben ser data.frames
  
  nr=nrow(sr)
  nv=nrow(sv)
  y=sv[, name_y]
  di=1/pi
  Nr=sum(di)
  if(pesos=="valliant")
    wi_ipw=valliant_weights(ps_sv) #CAMBIAR
  if(pesos=="sc")
    wi_ipw=sc_weights(ps_sv) #CAMBIAR
  if(pesos=="kernel")
    wi_ipw=kw.wt(p_score.c =ps_sv, p_score.s =ps_sr, svy.wt = N/n_A)$pswt #CAMBIAR 
  
  Nv=sum(wi_ipw)
  
  
  if(is.null(pij_mat)==TRUE){
    pij_mat=Pkl.Hajek.s(pi)
  }
  
  e=y-y_pred_sv;
  sigma2=var(e)*((length(e)-1)/length(e))
  
  ## Var1 ##
  #Va
  sum_sv=0
  sum_sr=0
  
  for(i in 1:nv){
    sum_sv=sum_sv+((1-2*ps_sv[i])/ps_sv[i]^2)*sigma2
  }
  #sum_sv
  
  for(i in 1:nr){
    sum_sr=sum_sr + di[i]*sigma2
  }
  va=sum_sv/(Nv^2) + sum_sr/(Nr^2)
  
  
  #Vb
  vb=0
  sum_vb=matrix(NA, nr, nr)
  for(i in 1:nr){
    sum_vb[i,]=parSapply(cluster, 1:nr,   sumatorio_mi_dr, i, pij_mat, di, y_pred_sr) #matriz nr*nr
  }
  vb=sum(sum_vb) 
  
  
  v_dr_b=vb/Nr^2
  v_dr1=va+v_dr_b
  
  
  ## Var2 ##
  #wa
  sum_sv=0
  for(i in 1:nv){
    sum_sv=sum_sv+((1-ps_sv[i])/ps_sv[i]^2)*(e[i]^2)
  }
  wa=sum_sv/Nv^2
  
  
  #wb
  wb=0
  sum_wb=matrix(NA, nr, nr)
  for(i in 1:nr){
    sum_wb[i,]=parSapply(cluster, 1:nr,   sumatorio_mi_dr, i, pij_mat, di, y_pred_sr) #matriz nr*nr
  }
  wb=sum(sum_wb) 
  v_wb=wb/Nr^2
  
  
  #wc
  sum_sv=0
  for(i in 1:nv){
    sum_sv=sum_sv+(sigma2/ps_sv[i])
  }
  
  sum_sr=0
  for(i in 1:nr){
    sum_sr=sum_sr+(di[i]*sigma2) #di tiene dimension nr 100, no nv 120---MIRAR
  }
  wc=sum_sv/(Nv^2)-sum_sr/(Nr^2)
  
  v_dr2=wa+v_wb-wc
  
  return(list(v1=v_dr1, v2=v_dr2))
}

