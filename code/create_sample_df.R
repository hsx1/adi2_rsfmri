library(dplyr) # version 1.0.2
library(tidyr) # version 1.1.2


create_sample_df <- function(group = "all", tp = "all"){
  # possible values for group: IG/KG/both
  # possible values for tp: BL/FU/FU2/BLFU/all
  
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
  full_sample=full_sample[,c("subj.ID","condition","tp","Age_BL","Sex","meanFD", "BMI","Exclude")]
  
  # Exclusion of participants --------------------------------------------------
  # exclusion due to problems in preprocessing
  full_sample$exclude_prep=full_sample$Exclude
  full_sample$exclude_prep[is.na(full_sample$exclude_prep)]=FALSE
  mri_sample=merge(mri_files, full_sample[!is.na(full_sample$condition),c("subj.ID","condition","tp","Age_BL","Sex","meanFD", "BMI","exclude_prep")], by=c("subj.ID","tp"))
  
  mri_sample[mri_sample$exclude_prep==TRUE,] # ADI063_bl ADI063_fu ADI116_bl ADI116_fu ADI116_fu2
  cat("Note that",nrow(mri_sample[mri_sample$exclude_prep==TRUE,]),"data points were excluded due to problems in preprocessing.\n")
  sample_after_1stExcl=mri_sample[mri_sample$exclude_prep==FALSE,]
  
  # exclusion if mean FD > average meanFD + 2 SD
  sample_after_1stExcl$exclude_fd = (sample_after_1stExcl$meanFD > mean(sample_after_1stExcl$meanFD)+ 2*sd(sample_after_1stExcl$meanFD))
  cat("Note that",nrow(sample_after_1stExcl[sample_after_1stExcl$exclude_fd==TRUE,]),"data points were excluded due to extensive head motion.\n")
  
  final_sample=sample_after_1stExcl[sample_after_1stExcl$exclude_fd==FALSE,]
  cat("The remaining sample comprises",nrow(final_sample),"data points across across time for all groups.\n")
  
  df=suppressMessages(plyr::join(mri_files,final_sample))
  rm(mri_files, full_sample, final_sample, mri_sample, tmp, i)
  
  
  # preparation of variables ---------------------------------------------------
  
  df$logmFD <- log10(df$meanFD)
  
  df$group[df$condition == "IG"] = 1
  df$group[df$condition == "KG"] = 2
  df$group_factor=as.factor(df$group)
  levels(df$group_factor)=c("IG","KG")
  
  df$tp=as.factor(df$tp)
  levels(df$tp)=c("bl","fu","fu2")
  
  df$visit=df$tp
  levels(df$visit)=c(1,2,3)
  
  df$Sex=as.factor(df$Sex) 
  levels(df$Sex)=c(-1,1)
  
  
  # Modify Design matrix for subsamples ----------------------------------------
  
  # Model specification incl. design matrix for different options
  # selection of groups
  if (group == "both") {
    # total: both groups = original dataframe
  }else if (group == "IG") {
    df <- df[df$condition=="IG",]
  }else if (group == "KG") {
    df <- df[df$condition=="KG",]
  }
  
  # selection of time points
  if (tp == "all"){
    df <- df
  }else if (tp == "BL"){
    df <- df[df$tp == "bl",]
  }else if (tp == "FU"){
    df <- df[df$tp == "fu",]
  }else if (tp == "FU2"){
    df <- df[df$tp == "fu2",]
  }else if (tp == "BLFU"){
    df <- df[df$tp != "fu2",]
  }else if (tp == "FUFU2"){
    df <- df[df$tp != "fu2",]
  }

  # ----------------------------------------------------------------------------
  # Additional computations (also better to implement in Matlab)
  # ----------------------------------------------------------------------------
  
  if ((group == "all" | group == "IG") & (tp == "all" | tp == "BLFU")) {
    # compute centered avgBMI (avgFD) and cgnBMI(cgnFD)
    df_wide <- tidyr::pivot_wider(
      data = df,
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
      # CAREFUL: Order according to dataframe df
      names_to = "tp_tp",
      values_to = "tp",
      values_drop_na = TRUE
    )
    all(df_long$tp == df$tp) # check if order correct
    
    head(df[, c("subj.ID", "tp")])
    head(df_long[, c("subj.ID", "tp")])
    
    df_long$cgnBMI <- df$BMI - df_long$avgBMI
    df_long$avgBMIc <- df_long$avgBMI - mean(df_long$avgBMI, na.rm = TRUE)
    df$avgBMIc <- df_long$avgBMIc
    df$cgnBMI <- df_long$cgnBMI
    
    df_long$cgnFD <- df$logmFD - df_long$avgFD
    df_long$avgFDc <- df_long$avgFD - mean(df_long$avgFD, na.rm = TRUE)
    df$avgFDc <- df_long$avgFDc
    df$cgnFD <- df_long$cgnFD
  }
  # end of extra computations --------------------------------------------------
  
  return(df)
}
