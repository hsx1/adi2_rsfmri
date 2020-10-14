library(dplyr)
library(tidyr)

save_txt_for_swe <- function(IG_only, only2tp){
  parentdir <- "/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /SwE_files"
  setwd(parentdir)
  
  files=read.table("scans_PCC_CC_z.txt")
  #split info from path in txt
  for (i in 1:nrow(files)){
    #print(strsplit(toString(rs_QA[i,1]),'_')[[1]][1])
    tmp=strsplit(toString(files[i,1]),'/')[[1]][8]
    files[i,"scan_dir"]=paste0(strsplit(toString(files[i,1]),'/')[[1]][1:8], sep = "/", collapse="")
    files[i,"tp"]=strsplit(toString(tmp),'_')[[1]][2]
    files[i,"subj.ID"]=strsplit(toString(tmp),'_')[[1]][1]
    files[i,"subj.Nr"]=as.numeric(substr(files$subj.ID[i],4,6))
  }
  files$V1 <- NULL
  
  #load group info and merge with path info
  info_file=read.csv("/data/p_02161/ADI_studie/metadata/final_sample_MRI_QA_info.csv")
  group_info=info_file[,c("subj.ID","condition","tp","Age_BL","Sex","meanFD", "BMI")]
  condition=merge(files, group_info[!is.na(group_info$condition),], by=c("subj.ID","tp"))
  
  final=plyr::join(files,condition)
  rm(files, group_info, info_file, condition, tmp, i)
  

# preparation of covariates -----------------------------------------------
  
  #coding "IG" "KG"
  final[final$condition=="IG","IG"]=1
  final[is.na(final$IG),"IG"]=0
  final[final$condition=="KG","KG"]=1
  final[is.na(final$KG),"KG"]=0
  
  # preparation of variables
  final$group[final$condition == "IG"] = 1
  final$group[final$condition == "KG"] = 2
  final$logmFD <- log10(final$meanFD)
  
  final$tp=as.factor(final$tp)
  final$visits=final$tp
  levels(final$visits)=c(1,2,3)
  final$tp_cov=final$tp
  levels(final$tp_cov)=c(-1,0,1)
  
  final$Sex=as.factor(final$Sex)
  levels(final$Sex)=c(-1,1)
  
  # linear modeling of time for each group
  final$tp_IG <- final$tp_cov
  final[final$KG==1,"tp_IG"]=0
  final$tp_KG <- final$tp_cov
  final[final$IG==1,"tp_KG"]=0
  
  # model with time as factor
  final[final$tp=="fu2"&final$IG==1,"IG_fu2"]=1
  final[final$tp=="fu"&final$IG==1,"IG_fu"]=1
  final[final$tp=="bl"&final$IG==1,"IG_bl"]=1
  final[is.na(final$IG_fu2),"IG_fu2"]=0
  final[is.na(final$IG_fu),"IG_fu"]=0
  final[is.na(final$IG_bl),"IG_bl"]=0
  
    final[final$tp=="fu2"&final$KG==1,"KG_fu2"]=1
    final[final$tp=="fu"&final$KG==1,"KG_fu"]=1
    final[final$tp=="bl"&final$KG==1,"KG_bl"]=1
    final[is.na(final$KG_fu2),"KG_fu2"]=0
    final[is.na(final$KG_fu),"KG_fu"]=0
    final[is.na(final$KG_bl),"KG_bl"]=0
  
  levels(final$tp)=c("bl","fu","fu2")

# Modify Design matrix for subsamples -------------------------------------

  # Model specification incl. design matrix for different options
  if (only2tp == TRUE){
    # only BL and FU - both groups
    output_dir <- file.path(parentdir, "only2tp")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    final <- final[final$tp != "fu2",]
    if (IG_only == TRUE){
      # only BL and FU - only IG
      if (!dir.exists(output_dir)) {dir.create(output_dir)}
      output_dir <- file.path(parentdir, "only2tp/IG_only")
      setwd(output_dir)
      final <- final[final$condition=="IG",]
    }
  }else if (only2tp == FALSE){
    # all time points
    if (IG_only == TRUE){
      # only IG
      output_dir <- file.path(parentdir, "IG_only")
      if (!dir.exists(output_dir)) {dir.create(output_dir)}
      setwd(output_dir)
      final <- final[final$condition=="IG",]
    }else{
      # total: both groups
      output_dir <- file.path(parentdir, "total")
      if (!dir.exists(output_dir)) {dir.create(output_dir)}
      setwd(output_dir)
    }
  }
  
  # directory load scans
  final[final$tp =="bl","scan_dir"]
  nrow(final[final$tp =="fu",])
  nrow(final[final$tp =="fu2",])
  write.table(final$scan_dir, col.names=FALSE,row.names=FALSE,quote=FALSE,
              file='scans.txt')
  # Modified SwE type - Visits: tp.txt
  write.table(final$visits, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='tp.txt')
  # Modified SwE type - Groups: group.txt
  write.table(final$group, col.names=FALSE,row.names=FALSE,quote=FALSE,
              file='group.txt')
  # Subjects
  write.table(final$subj.ID, col.names=FALSE, row.names=FALSE,
              file='subjID.txt')
  write.table(final$subj.Nr, col.names=FALSE, row.names=FALSE,
              file='subjNr.txt')
  
  ## Covariates
  # nuisance covariates
  write.table(final$Age_BL, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='Age.txt')
  write.table(final$Sex, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='Sex.txt')
  write.table(final$logmFD, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='logmeanFD.txt')
  # for covariates of interest
  write.table(final$IG, col.names=FALSE,row.names=FALSE, quote=FALSE,
              file='group_IG.txt')
  write.table(final$KG, col.names=FALSE,row.names=FALSE,quote=FALSE,
              file='group_KG.txt')
  
  #######################################
  # Different model with time as factor #
  #######################################
  
  write.table(final$tp_IG, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='tp_IG.txt')
  write.table(final$tp_KG, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='tp_KG.txt')
  write.table(final$IG_fu2, col.names=FALSE,row.names=FALSE,quote=FALSE,
              file='IG_fu2.txt')
  write.table(final$IG_fu, col.names=FALSE,row.names=FALSE,quote=FALSE,
              file='IG_fu.txt')
  write.table(final$IG_bl, col.names=FALSE,row.names=FALSE,quote=FALSE,
              file='IG_bl.txt')

    write.table(final$KG_fu2, col.names=FALSE,row.names=FALSE,quote=FALSE,
              file='KG_fu2.txt')
    write.table(final$KG_fu, col.names=FALSE,row.names=FALSE,quote=FALSE,
              file='KG_fu.txt')
    write.table(final$KG_bl, col.names=FALSE,row.names=FALSE,quote=FALSE,
              file='KG_bl.txt')

  
  # compute centered avgBMI (avgFD) and cgnBMI(cgnFD)
  df_wide <- tidyr::pivot_wider(data = final,
                                id_cols = "subj.ID",
                                names_from = "tp",
                                values_from = c("tp","BMI","logmFD")
  )
  
  if (only2tp == TRUE) {s = 2} else {s = 1}
  
  BMI_vector <- c("BMI_fu2", "BMI_fu", "BMI_bl")
  FD_vector <- c("logmFD_fu2", "logmFD_fu", "logmFD_bl")
  df_wide$avgBMI <- apply(df_wide[, BMI_vector[s:3]],
                          MARGIN = 1,
                          FUN = mean,
                          na.rm = TRUE)
  df_wide$avgFD <- apply(df_wide[, FD_vector[s:3]],
                         MARGIN = 1,
                         FUN = mean,
                         na.rm = TRUE)
  
  if (only2tp == TRUE) {order_vector = c(3,2)} else {order_vector = c(4, 2, 3)}
  df_long <- tidyr::pivot_longer(
    data = df_wide,
    cols = all_of(order_vector), # bl, fu, fu2 (maintain original order)
    # CAREFUL: Order according to dataframe final
    names_to = "tp_tp",
    values_to = "tp",
    values_drop_na = TRUE
  )
  all(df_long$tp == final$tp) # check if order correct
  
  head(final[,c("subj.ID","tp")])
  head(df_long[,c("subj.ID","tp")])
  
  
  df_long$cgnBMI <- final$BMI - df_long$avgBMI
  df_long$avgBMIc <- df_long$avgBMI - mean(df_long$avgBMI,na.rm=TRUE)
  final$avgBMIc <- df_long$avgBMIc
  final$cgnBMI <- df_long$cgnBMI
  
  df_long$cgnFD <- final$logmFD - df_long$avgFD
  df_long$avgFDc <- df_long$avgFD - mean(df_long$avgFD,na.rm=TRUE)
  final$avgFDc <- df_long$avgFDc
  final$cgnFD <- df_long$cgnFD
  
  write.table(final$avgBMIc, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='avgBMIc.txt')
  write.table(final$cgnBMI, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='cgnBMI.txt')
  write.table(final$BMI, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='BMI.txt')
  write.table(final$avgFDc, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='avgFDc.txt')
  write.table(final$cgnBMI, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='cgnFD.txt')
  write.table(final$meanFD, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='meanFD.txt')
  return(final)
}



# then save txt for all groups and return the final dataframe
save_txt_for_swe(IG_only = FALSE, only2tp = FALSE)
# save txt first only IG (under different directory as specified in the function)
# save_txt_for_swe(IG_only = TRUE)


