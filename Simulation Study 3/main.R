setwd("...")

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

cluster = makeCluster(3) #para que no coja todos los nucleos? OJO QUE LOS COJO TODOS PARA QUE TARDE MENOS
clusterEvalQ(cluster,library(MASS))
clusterEvalQ(cluster,library(tidyverse))
clusterEvalQ(cluster,library(mlbench))
clusterEvalQ(cluster,library(xgboost))
clusterEvalQ(cluster,library(foreign))
clusterEvalQ(cluster,library(RcmdrMisc))
clusterEvalQ(cluster,library(fastDummies))
clusterEvalQ(cluster,library(lme4))
clusterEvalQ(cluster,library(NonProbEst))
clusterEvalQ(cluster,library(caret))
clusterEvalQ(cluster,library(sampling))
clusterEvalQ(cluster,library(survey))
#clusterEvalQ(cluster,library(FSelector))
clusterEvalQ(cluster,library(KWML))
clusterEvalQ(cluster,library(samplingVarEst))
#remotes::install_github("ncn-foreigners/nonprobsvy")
clusterEvalQ(cluster,library(nonprobsvy))
clusterEvalQ(cluster,library(dplyr))
clusterEvalQ(cluster,library(Frames2))
#clusterEvalQ(cluster,library(xlsx))
clusterEvalQ(cluster,library(KWML))
clusterEvalQ(cluster,library(nnet))


## Variances 
# Feliz
source("simulation feliz.R") 
setwd("...") 
write.xlsx(results_feliz, "results_feliz.xlsx", sheetName="GLM")

# Percepción
source("simulation perc.R") 
setwd("...")
write.xlsx(results_perc, "results_perc.xlsx", sheetName="GLM")


## ML scenario
# Feliz 
source("simulation feliz ML.R") 
setwd("...")
write.xlsx(results_feliz, "results_feliz_ML.xlsx", sheetName="GLM")
write.xlsx(knn_feliz, "results_feliz_ML.xlsx", sheetName="KNN", append=TRUE)
write.xlsx(gbm_feliz, "results_feliz_ML.xlsx", sheetName="GBM", append=TRUE)
write.xlsx(crf_feliz, "results_feliz_ML.xlsx", sheetName="CRF", append=TRUE)

# Percepción
source("simulation perc ML.R") 
setwd("...") 
write.xlsx(results_perc, "results_perc_ML.xlsx", sheetName="GLM")
write.xlsx(knn_perc, "results_perc_ML.xlsx", sheetName="KNN", append=TRUE)
write.xlsx(gbm_perc, "results_perc_ML.xlsx", sheetName="GBM", append=TRUE)
write.xlsx(crf_perc, "results_perc_ML.xlsx", sheetName="CRF", append=TRUE)


## Scenario para nr>nv
# Feliz 
source("simulation feliz nr.R") 
setwd("...")
write.xlsx(results_feliz, "results feliz nr.xlsx", sheetName="GLM")

# Percepción
source("simulation perc nr.R") 
setwd("...")
write.xlsx(results_perc, "results perc nr.xlsx", sheetName="GLM")

stopCluster()