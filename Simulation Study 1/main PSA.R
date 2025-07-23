setwd("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Escenario Inicial (estudio simul art)/nr=1000 PSA")

#for(n_r in c(500, 1000, 1500, 2000, 2500)){
  n_r=1000
  source("simulation PSA1000.R")  
  setwd("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Escenario Inicial (estudio simul art)/nr=1000 PSA")
  write.xlsx(results, paste0("results_PSA1000_",n_r,".xlsx"))
#}



# setwd("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Escenario Inicial (estudio simul art)/datos")
# nv=vector()
# for(i in 1:rep){
#   setwd(file.path(paste0("Sample ",i)))
#   sv=read.table(file="muestra50repl.txt", header = FALSE, sep = "", dec = ".") #muestra non-prob (RDS)
#   nv[i]=nrow(sv)
#   setwd("C:/Users/jorge/Desktop/EST/TI/2024-2025 (2º Año Doctorado)/Artículos/Artículo Ismael/Escenario Inicial (estudio simul art)/datos")
# }
# 
# summary(nv)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 215.0   238.8   251.0   250.3   261.0   285.0 
