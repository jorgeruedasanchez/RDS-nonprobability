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
library(mice)

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

meanest <- function(variable, s_v, wi_i) svymean(as.formula(paste0("~", variable)), 
                                                 design = trimWeights(svydesign(~0, data = s_v, weights = wi_i), lower = 1/sum(wi_i)), 
                                                 method = "mean")

matching_prop <- function (convenience_sample, reference_sample, covariates, estimated_var, 
                      positive_label = NULL, algorithm, proc = NULL,
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

matching_num <- function (convenience_sample, reference_sample, covariates, estimated_var, 
                          algorithm, proc = NULL,
                          trControl = trainControl(method = "none"), ...) 
{
  data = convenience_sample[, covariates, drop = FALSE]
  values = convenience_sample[, estimated_var]
  test = reference_sample[, covariates, drop = FALSE]
  model = train(data, values, algorithm, preProcess = proc, 
                trControl = trControl, ...)
    return(predict(model, test))
}

matching_num_xgb <- function (convenience_sample, reference_sample, covariates, estimated_var, 
                      positive_label = NULL, algorithm = "xgbTree", proc = NULL,
                      trControl = trainControl(method = "none"), ...) 
{
  data = convenience_sample[, covariates, drop = FALSE]
  values = convenience_sample[, estimated_var]
  test = reference_sample[, covariates, drop = FALSE]
  model = train(data, values, algorithm, preProcess = proc, trControl = trControl, 
                tuneGrid = data.frame(nrounds=100, max_depth=6, eta=0.1, gamma=0, colsample_bytree=1, min_child_weight=1, subsample=1),...)
    return(predict(model, test))
}


model_todos <- function(sample_data, weights, full_data, covariates, estimated_var, 
                        estimate_mean = FALSE, positive_label = NULL, algorithm, 
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

propensities <- function (convenience_sample, reference_sample, covariates, algorithm, model_weights = NULL,
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


propensities_xgb <- function (convenience_sample, reference_sample, covariates, algorithm = "xgbTree", model_weights = NULL,
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
                                                                  "Negative")), algorithm, weights = model_weights, preProcess = proc, trControl = trControl, 
                tuneGrid = data.frame(nrounds=100, max_depth=6, eta=0.1, gamma=0, colsample_bytree=1, min_child_weight=1, subsample=1),...)
  probabilities = predict(model, data, type = "prob")$Positive
  if (smooth) 
    probabilities = (1000 * probabilities + 0.5)/1001
  list(convenience = probabilities[1:n_convenience], reference = probabilities[(n_convenience + 
                                                                                  1):length(probabilities)])
}

papp <- function (convenience_sample, reference_sample, covariates, weights, 
                  algorithm, proc = NULL, ...) 
{
  test = convenience_sample[, covariates, drop = FALSE]
  values = weights
  data = reference_sample[, covariates, drop = FALSE]
  model = train(data, values, algorithm, preProcess = proc, 
                trControl = trainControl(classProbs = TRUE, method = "none"), ...)
  return(predict(model, test))
}

psaplusmatch <- function(muestra_cohorte, muestra_prob, pesos, di, X, algoritmo,
                             trControl = trainControl(classProbs = TRUE, method = "none")){
  npesos <- nrow(muestra_cohorte)*pesos/sum(pesos)
  y_matching_psa <- matching(muestra_cohorte %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))),
                             muestra_prob %>% mutate(y = factor(y, levels = c(0, 1), labels = c("X0", "X1"))), 
                             X, "y", positive_label = "X1", algorithm = algoritmo, weights = npesos,
                             trControl = trControl)
  return(c(sum(y_matching_psa * di)/3000000, sum(y_matching_psa * di)/sum(di)))
}

