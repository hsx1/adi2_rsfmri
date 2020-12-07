library(dplyr) # version 1.0.2
library(tidyr) # version 1.1.2

create_sample_df <- function(group = "all", tp = "all"){
  # possible values for group: IG/KG/both
  # possible values for tp: BL/FU/FU2/BLFU/all
  
  # participants with mri data (txt file with path to .nii file) ---------------
  
  mri_files=read.table("../SwE_files/scans_PCC_CC_z.txt")
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
  
  
  # preparation of variables ---------------------------------------------------
  
  final$logmFD <- log10(final$meanFD)
  
  final$group[final$condition == "IG"] = 1
  final$group[final$condition == "KG"] = 2
  final$group_factor=as.factor(final$group)
  levels(final$group_factor)=c("IG","KG")
  
  final$tp=as.factor(final$tp)
  levels(final$tp)=c("bl","fu","fu2")
  
  final$visit=final$tp
  levels(final$visit)=c(1,2,3)
  
  final$Sex=as.factor(final$Sex) 
  levels(final$Sex)=c(-1,1)
  
  
  # Modify Design matrix for subsamples ----------------------------------------
  
  # Model specification incl. design matrix for different options
  # selection of groups
  if (group == "both") {
    # total: both groups = original dataframe
  }else if (group == "IG") {
    final <- final[final$condition=="IG",]
  }else if (group == "KG") {
    final <- final[final$condition=="KG",]
  }
  
  # selection of time points
  if (tp == "all"){
    final <- final
  }else if (tp == "BL"){
    final <- final[final$tp == "bl",]
  }else if (tp == "FU"){
    final <- final[final$tp == "fu",]
  }else if (tp == "FU2"){
    final <- final[final$tp == "fu2",]
  }else if (tp == "BLFU"){
    final <- final[final$tp != "fu2",]
  }else if (tp == "FUFU2"){
    final <- final[final$tp != "fu2",]
  }

  # ----------------------------------------------------------------------------
  # Additional computations (also better to implement in Matlab)
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
  }
  # end of extra computations --------------------------------------------------
  
  return(final)
}
