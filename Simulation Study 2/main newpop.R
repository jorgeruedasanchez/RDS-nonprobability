library(MASS)
library(tidyverse)
library(mlbench)
library(caret)
library(survey)
library(sampling)
library(NonProbEst)
library(xgboost)
library(dplyr)
library(foreign)
library(RcmdrMisc)
library(samplingVarEst)
# install.packages("devtools")
# library(devtools)
library(fastDummies)
library(parallel)
library(lme4)
#library(FSelector)
#remotes::install_github("ncn-foreigners/nonprobsvy")
library(nonprobsvy)
library(Frames2)
#library(xlsx)
# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("chkern/KWML")
library(KWML)
#install.packages("nnet")
library(nnet)
library(xlsx)

n_redes1=c("Highly Connected Groups", "Scale free", "Small Groups")
n_r=1000 #Como en el que calculamos las varianzas

for(escenario in n_redes1){
  source("simulation newpop1.R")  
  setwd("...")
  write.xlsx(results, paste0("results",escenario,".xlsx"))
}


n_redes2=c("sim2a 3 seeds", "sim2b 8 seeds")

for(escenario in n_redes2){
  source("simulation newpop2.R")  
  setwd("...")
  write.xlsx(results, paste0("results",escenario,".xlsx"))
}
