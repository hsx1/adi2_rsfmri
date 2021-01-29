library(dplyr) # version 1.0.2
library(tidyr) # version 1.1.2


create_sample_df <- function(group = "both", tp = "all", exclFD = FALSE){
  # possible values for group: IG/KG/both
  # possible values for tp: BL/FU/FU2/BLFU/all

  # ----------------------------------------------------------------------------
  # second level functions
  # ----------------------------------------------------------------------------

  apply_excl_criteria <- function(data = full_sample, exclFD = FALSE){
    # exclusion of participants who failed preprocessing in FreeSurfer
    cat("Fill sample consistst of N =",nrow(full_sample),"; ")
    full_sample$exclude_prep <- full_sample$Exclude
    full_sample$exclude_prep[is.na(full_sample$exclude_prep)]=FALSE
    
    mri_sample <- merge(mri_files, full_sample[!is.na(full_sample$condition),c("subj.ID","condition","tp","Age_BL","Sex","meanFD", "BMI","exclude_prep")], by=c("subj.ID","tp"))
    cat("n =",nrow(mri_sample),"with fMRI data.\n")

    mri_sample[mri_sample$exclude_prep==TRUE,] # ADI063_bl ADI063_fu ADI116_bl ADI116_fu ADI116_fu2
    cat("Note that",nrow(mri_sample[mri_sample$exclude_prep==TRUE,]),"data points were excluded due to problems in preprocessing.\n")
    sample_after_1stExcl <- mri_sample[mri_sample$exclude_prep==FALSE,]
    sample_after_1stExcl$exclude_prep <- NULL

    if (exclFD == TRUE){
      # exclusion of 10% of the worst mean FD
      thres10=sort(sample_after_1stExcl$meanFD, decreasing = TRUE)[round(nrow(sample_after_1stExcl)/10)]
      sample_after_1stExcl$exclude_fd=(sample_after_1stExcl$meanFD>=thres10)
      cat("Note that",nrow(sample_after_1stExcl[sample_after_1stExcl$exclude_fd == TRUE, ]),"data points were excluded due to extensive head motion.\n")

      final_sample=sample_after_1stExcl[sample_after_1stExcl$exclude_fd==FALSE,]
      sample_after_1stExcl$exclude_fd <- NULL
    }else if (exclFD ==FALSE){
      final_sample=sample_after_1stExcl
    }
    cat("The remaining sample comprises n =",nrow(final_sample),"data points across across time for both groups.\n")
    # all remaining data points should be included
    final_sample$include <- TRUE
    return(final_sample)
  }

  # ----------------------------------------------------------------------------
  # script
  # ----------------------------------------------------------------------------

  # participants with mri data (txt file with path to .nii file) ---------------

  mri_files=data.frame(list.files("/data/pt_02161/Results/Project2_resting_state/connectivity/calc_DMN_reward_seed_connectivity",
             full.names = TRUE))
  
  # split info from path in txt
  for (i in 1:nrow(mri_files)){
    tmp=strsplit(toString(mri_files[i,1]),'/')[[1]][8]
    mri_files[i,"scan_dir"]=paste0(strsplit(toString(mri_files[i,1]),'/')[[1]][1:8], sep = "/", collapse="")
    mri_files[i,"tp"]=strsplit(toString(tmp),'_')[[1]][2]
    mri_files[i,"subj.ID"]=strsplit(toString(tmp),'_')[[1]][1]
    mri_files[i,"subj.Nr"]=as.numeric(substr(mri_files$subj.ID[i],4,6))
  }
  mri_files[,1] <- NULL

  # load info from full sample and merge with path info (fmri only)-------------

  full_sample=read.csv("/data/p_02161/ADI_studie/metadata/final_sample_MRI_QA_info.csv") # CAVE: info file elsewhere !!!
  full_sample=full_sample[,c("subj.ID","condition","tp","Age_BL","Sex","meanFD","maxFD", "BMI", "Final_Score",
                           "Excessive_motion", "Exclude", "Check_comments")]

  # Exclusion of participants --------------------------------------------------
  final_sample <- apply_excl_criteria(full_sample, exclFD)

  df=suppressMessages(plyr::join(mri_files,final_sample))
  df=df[!is.na(df$include),]
  rm(apply_excl_criteria, mri_files, full_sample, final_sample, tmp, i)

  # preparation of variables --------------------------

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
  cat("Data from",group,"group(s) and",tp,"measurements were selected to be kept in the data frame, therefore n =",nrow(df),".\n")

  # end of computations --------------------------------------------------------
  return(df)
}
