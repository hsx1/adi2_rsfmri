library(dplyr) # version 1.0.2
library(tidyr) # version 1.1.2

get_txt_for_swe <- function(group = "all", tp = "all"){
  # group = IG/KG/both
  # tp = BL/FU/FU2/BLFU/FUFU2/all
  
  original_wd=getwd()
  # import absolute path
  abs_path=read.csv("abs_path.csv", header=FALSE, stringsAsFactors=FALSE)
  # set working directory
  parentdir=file.path(getwd(), "../SwE_files/", fsep = .Platform$file.sep)[1]
  setwd(parentdir)
  
  # participants with mri data (txt file with path to .nii file) ---------------
  
  mri_files=read.table("../SwE_files/scans.txt")
  # split info from path in txt
  for (i in 1:nrow(mri_files)){
    tmp=strsplit(toString(mri_files[i,1]),'/')[[1]][8]
    mri_files[i,"scan_dir"]=paste0(strsplit(toString(mri_files[i,1]),'/')[[1]][1:8], sep = "/", collapse="")
    mri_files[i,"tp"]=strsplit(toString(tmp),'_')[[1]][2]
    mri_files[i,"subj.ID"]=strsplit(toString(tmp),'_')[[1]][1]
    mri_files[i,"subj.Nr"]=as.numeric(substr(mri_files$subj.ID[i],4,6))
  }
  mri_files$V1 <- NULL
  
  
  # load info from full sample and merge with path info (fmri only)-------------
  
  full_sample=read.csv("/data/p_02161/ADI_studie/metadata/final_sample_MRI_QA_info.csv") # CAVE: info file elsewhere !!!
  full_sample=full_sample[,c("subj.ID","condition","tp","Age_BL","Sex","meanFD", "BMI")]
  condition=merge(mri_files, full_sample[!is.na(full_sample$condition),], by=c("subj.ID","tp"))
  final=plyr::join(mri_files,condition)
  rm(mri_files, full_sample, condition, tmp, i)
  
# preparation of covariates ----------------------------------------------------
  
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
  final$visit=final$tp
  levels(final$visit)=c(1,2,3)
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
  output_dir <- parentdir
  # selection of groups
  if (group == "both") {
    # total: both groups
  }else if (group == "IG") {
    output_dir <- file.path(parentdir, "IG_only")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    final <- final[final$condition=="IG",]
  }else if (group == "KG") {
    output_dir <- file.path(parentdir, "KG_only")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    final <- final[final$condition=="KG",]
  }
  
  # selection of time points
  if (tp == "all"){
    output_dir <- file.path(output_dir, "total")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    final <- final
  }else if (tp == "BL"){
    output_dir <- file.path(output_dir, "onlyBL")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    final <- final[final$tp == "bl",]
  }else if (tp == "FU"){
    output_dir <- file.path(output_dir, "onlyFU")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    final <- final[final$tp == "fu",]
  }else if (tp == "FU2"){
    output_dir <- file.path(output_dir, "onlyFU2")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    final <- final[final$tp == "fu2",]
  }else if (tp == "BLFU"){
    output_dir <- file.path(output_dir, "only2tp")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    final <- final[final$tp != "fu2",]
  }else if (tp == "FUFU2"){
    output_dir <- file.path(output_dir, "only2tp_fu")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    final <- final[final$tp != "fu2",]
  }

  # directory load scans
  write.table(final$scan_dir, col.names=FALSE,row.names=FALSE,quote=FALSE,
              file='scans.txt')
  # Modified SwE type - visit: tp.txt
  write.table(final$visit, col.names=FALSE, row.names=FALSE,quote=FALSE,
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
  
  # other important variables
  write.table(final$BMI, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='BMI.txt')
  write.table(final$meanFD, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='meanFD.txt')
  
  # ----------------------------------------------------------------------------
  # Different model with time as factor (instead of continous)
  # ----------------------------------------------------------------------------
  
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

  # ----------------------------------------------------------------------------
  # Extra analysis steps (also better implmented in Matlab)
  # ----------------------------------------------------------------------------
    
  if ((group == "all" | group == "IG") & (tp == "all" | tp == "BLFU")) {
      # compute centered avgBMI (avgFD) and cgnBMI(cgnFD)
      df_wide <- tidyr::pivot_wider(
        data = final,
        id_cols = "subj.ID",
        names_from = "tp",
        values_from = c("tp", "BMI", "logmFD")
      )
      
      if (tp == "BLFU") {s = 2} else {s = 1}
      
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
      
      if (tp == "BLFU") {
        order_vector = c(3, 2)
      } else {
        order_vector = c(4, 2, 3)
      }
    
      df_long <- tidyr::pivot_longer(
        data = df_wide,
        cols = all_of(order_vector), # bl, fu, fu2 (maintain original order)
        # CAREFUL: Order according to dataframe final
        names_to = "tp_tp",
        values_to = "tp",
        values_drop_na = TRUE
      )
      all(df_long$tp == final$tp) # check if order correct
      
      head(final[, c("subj.ID", "tp")])
      head(df_long[, c("subj.ID", "tp")])
      
      df_long$cgnBMI <- final$BMI - df_long$avgBMI
      df_long$avgBMIc <- df_long$avgBMI - mean(df_long$avgBMI, na.rm = TRUE)
      final$avgBMIc <- df_long$avgBMIc
      final$cgnBMI <- df_long$cgnBMI
      
      df_long$cgnFD <- final$logmFD - df_long$avgFD
      df_long$avgFDc <- df_long$avgFD - mean(df_long$avgFD, na.rm = TRUE)
      final$avgFDc <- df_long$avgFDc
      final$cgnFD <- df_long$cgnFD
      
      write.table(final$avgBMIc, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='avgBMIc.txt')
      write.table(final$cgnBMI, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='cgnBMI.txt')
      write.table(final$avgFDc, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='avgFDc.txt')
      write.table(final$cgnBMI, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='cgnFD.txt')
  }

  # reset to original working directory
  setwd(original_wd)
  
}

# then save txt for all groups and return the final dataframe
get_txt_for_swe(group = "all", tp = "all")
