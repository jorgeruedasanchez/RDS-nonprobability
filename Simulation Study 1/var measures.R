 load("results_var.Rdata")
 RV=apply(est, 2, var) 
 SD=apply(est, 2, sd)

 #SE RDS
 esd1_RDSI=read.table("z_standard_error CI 1 RDS-I.txt")
 esd1_RDSII=read.table("z_standard_error CI 1 RDS-II.txt")
 esd1_RDSSS=read.table("z_standard_error CI 1 RDS-SS.txt")
 esd1_RDS=data.frame(esd1_RDSI, esd1_RDSII, esd1_RDSSS)
 colnames(esd1_RDS)=c("RDS-I", "RDS-II", "RDS-SS")
 SE1_RDS=apply(esd1_RDS, 2, mean)
 
 esd2_RDSI=read.table("z_standard_error CI 2 RDS-I.txt")
 esd2_RDSII=read.table("z_standard_error CI 2 RDS-II.txt")
 esd2_RDSSS=read.table("z_standard_error CI 2 RDS-SS.txt")
 esd2_RDS=data.frame(esd2_RDSI, esd2_RDSII, esd2_RDSSS)
 colnames(esd2_RDS)=c("RDS-I", "RDS-II", "RDS-SS")
 SE2_RDS=apply(esd2_RDS, 2, mean)
 
 #SE NP estimators
 esd_np=sqrt(dfvar)
 SE_np=apply(esd_np, 2, mean) 
 names(SE_np)=c("PSA1", "MI", "DR")
 SE1=c(SE1_RDS, SE_np)
 SE2=c(SE2_RDS, SE_np)

 RB1_sd=c(mean(abs(esd1_RDS[,1]-SD[1])/SD[1]),
         mean(abs(esd1_RDS[,2]-SD[2])/SD[2]),
         mean(abs(esd1_RDS[,3]-SD[3])/SD[3]),
         mean(abs(esd_np[,1]-SD[4])/SD[4]),
         mean(abs(esd_np[,2]-SD[5])/SD[5]),
         mean(abs(esd_np[,3]-SD[6])/SD[6]))*100 

 
 RRMSE1_sd=c(mean((esd1_RDS[,1]-SD[1])^2/SD[1]^2),
            mean((esd1_RDS[,2]-SD[2])^2/SD[2]^2),
            mean((esd1_RDS[,3]-SD[3])^2/SD[3]^2),
            mean((esd_np[,1]-SD[4])^2/SD[4]^2),
            mean((esd_np[,2]-SD[5])^2/SD[5]^2),
            mean((esd_np[,3]-SD[6])^2/SD[6]^2))^(1/2)
 
 RB2_sd=c(mean(abs(esd2_RDS[,1]-SD[1])/SD[1]),
          mean(abs(esd2_RDS[,2]-SD[2])/SD[2]),
          mean(abs(esd2_RDS[,3]-SD[3])/SD[3]),
          mean(abs(esd_np[,1]-SD[4])/SD[4]),
          mean(abs(esd_np[,2]-SD[5])/SD[5]),
          mean(abs(esd_np[,3]-SD[6])/SD[6]))*100 
 
 RRMSE2_sd=c(mean((esd2_RDS[,1]-SD[1])^2/SD[1]^2),
             mean((esd2_RDS[,2]-SD[2])^2/SD[2]^2),
             mean((esd2_RDS[,3]-SD[3])^2/SD[3]^2),
             mean((esd_np[,1]-SD[4])^2/SD[4]^2),
             mean((esd_np[,2]-SD[5])^2/SD[5]^2),
             mean((esd_np[,3]-SD[6])^2/SD[6]^2))^(1/2) 

 RBMSE1_df=data.frame(RB1_sd, RRMSE1_sd)
 rownames(RBMSE1_df)=c("RDS-I", "RDS-II", "RDS-SS", "PSA1", "MI", "DR")
 colnames(RBMSE1_df)=c("%RB", "RRMSE")
 
 SDSERB1_df=data.frame(SD, SE1, RB1_sd)
 rownames(SDSERB1_df)=c("RDS-I", "RDS-II", "RDS-SS", "PSA1", "MI", "DR")
 colnames(SDSERB1_df)=c("SD", "SE", "%RB")
 
 RBMSE2_df=data.frame(RB2_sd, RRMSE2_sd)
 rownames(RBMSE2_df)=c("RDS-I", "RDS-II", "RDS-SS", "PSA1", "MI", "DR")
 colnames(RBMSE2_df)=c("%RB", "RRMSE")
 
 SDSERB2_df=data.frame(SD, SE2, RB2_sd)
 rownames(SDSERB2_df)=c("RDS-I", "RDS-II", "RDS-SS", "PSA1", "MI", "DR")
 colnames(SDSERB2_df)=c("SD", "SE", "%RB")
  
 # install.packages("openxlsx")
 library(openxlsx)
 wb <- createWorkbook()
 addWorksheet(wb, "IC1 RDS")
 writeData(wb, sheet = "IC1 RDS", SDSERB1_df, startRow = 1, startCol = 1, rowNames = TRUE)
 writeData(wb, sheet = "IC1 RDS", RBMSE1_df, startRow = nrow(SDSERB1_df) + 3, startCol = 1, rowNames = TRUE)
 saveWorkbook(wb, "medvar.xlsx", overwrite = TRUE)
 
 addWorksheet(wb, "IC2 RDS")
 writeData(wb, sheet = "IC2 RDS", SDSERB2_df, startRow = 1, startCol = 1, rowNames = TRUE)
 writeData(wb, sheet = "IC2 RDS", RBMSE2_df, startRow = nrow(SDSERB2_df) + 3, startCol = 1, rowNames = TRUE)
 saveWorkbook(wb, "medvar.xlsx", overwrite = TRUE)