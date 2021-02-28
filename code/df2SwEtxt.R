library(dplyr) # version 1.0.2
library(tidyr) # version 1.1.2setwd()

get_txt_for_swe <- function(group = "both", tp = "all",exclFD = FALSE){
  # group = IG/KG/both
  # tp = BL/FU/FU2/BLFU/FUFU2/all
  
  # ----------------------------------------------------------------------------
  # primary definitions
  # ----------------------------------------------------------------------------
  
  original_wd=getwd()
  # import absolute path
  abs_path=read.csv("abs_path.csv", header=FALSE, stringsAsFactors=FALSE)
  # set working directory
  parentdir=file.path(getwd(), "../SwE_files", fsep = .Platform$file.sep)[1]
  setwd(parentdir)
  parentdir=getwd()
  
  
  # ----------------------------------------------------------------------------
  # second level functions
  # ----------------------------------------------------------------------------

  apply_excl_criteria <- function(data = full_sample, exclFD = FALSE){
    # exclusion due to problems in preprocessing
    cat("N=",nrow(full_sample),"\n")
    full_sample$exclude_prep <- full_sample$Exclude
    full_sample$exclude_prep[is.na(full_sample$exclude_prep)]=FALSE
    mri_sample <- merge(mri_files, full_sample[!is.na(full_sample$condition),c("subj.ID","condition","tp","Age_BL","Sex","meanFD", "BMI","exclude_prep")], by=c("subj.ID","tp"))
    cat("n=",nrow(mri_sample),"with fMRI data points.\n")
    
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
    cat("The remaining sample comprises",nrow(final_sample),"data points across across time for all groups.\n")
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
  
  # load info from full sample and merge with path info (fmri only) ------------
  
  full_sample=read.csv("/data/p_02161/ADI_studie/metadata/final_sample_MRI_QA_info.csv") # CAVE: info file elsewhere !!!
  full_sample=full_sample[,c("subj.ID","condition","tp","Age_BL","Sex","meanFD", "BMI","Exclude")]
  
  # Exclusion of participants --------------------------------------------------
  final_sample <- apply_excl_criteria(full_sample, exclFD)
  
  df=suppressMessages(plyr::join(mri_files,final_sample))
  df=df[!is.na(df$include),]
  rm(apply_excl_criteria, mri_files, full_sample, final_sample, tmp, i)

  
# preparation of covariates ----------------------------------------------------
  
  # coding "IG" "KG"
  df$IG[df$condition == "IG"] = 1
  df[is.na(df$IG), "IG"] = 0
  df$KG[df$condition == "KG"] = 1
  df[is.na(df$KG), "KG"] = 0
  
  # preparation of variables
  df$group[df$condition == "IG"] = 1
  df$group[df$condition == "KG"] = 2
  df$logmFD <- log10(df$meanFD)
  
  df$tp=as.factor(df$tp)
  df$visit=df$tp
  levels(df$visit)=c(1,2,3)
  df$tp_cov=df$tp
  levels(df$tp_cov)=c(-1,0,1)
  
  df$Sex=as.factor(df$Sex)
  levels(df$Sex)=c(-1,1)
  
  # linear modeling of time for each group
  df$tp_IG <- df$tp_cov
  df[df$KG==1,"tp_IG"]=0
  df$tp_KG <- df$tp_cov
  df[df$IG==1,"tp_KG"]=0
  
  # model with time as factor
  df[df$tp=="fu2"&df$IG==1,"IG_fu2"]=1
  df[df$tp=="fu"&df$IG==1,"IG_fu"]=1
  df[df$tp=="bl"&df$IG==1,"IG_bl"]=1
  df[is.na(df$IG_fu2),"IG_fu2"]=0
  df[is.na(df$IG_fu),"IG_fu"]=0
  df[is.na(df$IG_bl),"IG_bl"]=0
  
  df[df$tp=="fu2"&df$KG==1,"KG_fu2"]=1
  df[df$tp=="fu"&df$KG==1,"KG_fu"]=1
  df[df$tp=="bl"&df$KG==1,"KG_bl"]=1
  df[is.na(df$KG_fu2),"KG_fu2"]=0
  df[is.na(df$KG_fu),"KG_fu"]=0
  df[is.na(df$KG_bl),"KG_bl"]=0
  
  levels(df$tp)=c("bl","fu","fu2")

  df$BMIbl <- NA
  df$BMIbl[df$tp == "bl"] <- df$BMI[df$tp == "bl"]
  for (i in 1:nrow(df)){
    # check if bl exist for subject
    if (any(df$tp[df$subj.ID == df$subj.ID[i]] == "bl")){
      # replace with value of same subject at bl
      df$BMIbl[i] <- df$BMIbl[df$subj.ID == df$subj.ID[i] & df$tp =="bl"]
    }
  }
  
# Modify Design matrix for subsamples -------------------------------------
  
  # Model specification incl. design matrix for different options
  output_dir <- parentdir
  if (exclFD == FALSE){
    output_dir <- file.path(parentdir, "noExclFD") 
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
  }else if (exclFD == TRUE){
    output_dir <- file.path(parentdir, "ExclFD")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
  }
  
  # selection of groups
  if (group == "both"){
    output_dir <- file.path(output_dir, "both")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
  }else if (group == "IG"){
    output_dir <- file.path(output_dir, "IG")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    df <- df[df$condition=="IG",]
  }else if (group == "KG"){
    output_dir <- file.path(output_dir, "KG")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    df <- df[df$condition=="KG",]
  }
  
  # selection of time points
  
  if (tp == "all"){
    output_dir <- file.path(output_dir, "total")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    df <- df
  }else if (tp == "BL"){
    output_dir <- file.path(output_dir, "BL")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    df <- df[df$tp == "bl",]
  }else if (tp == "FU"){
    output_dir <- file.path(output_dir, "FU")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    df <- df[df$tp == "fu",]
  }else if (tp == "FU2"){
    output_dir <- file.path(output_dir, "FU2")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    df <- df[df$tp == "fu2",]
  }else if (tp == "BLFU"){
    output_dir <- file.path(output_dir, "BLFU")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    df <- df[df$tp != "fu2",]
  }else if (tp == "FUFU2"){
    output_dir <- file.path(output_dir, "FUFU2")
    if (!dir.exists(output_dir)) {dir.create(output_dir)}
    setwd(output_dir)
    df <- df[df$tp != "fu2",]
  }
  cat("Data from",group,"group(s) and",tp,"measurements were selected to be kept in the data frame, therefore n =",nrow(df),".\n")
  cat("All txt files will be saved in:", output_dir)
  
  
  # order dataframe by group / subject / time point
  df$tp_factor <- as.numeric(df$tp)
  dfo <- df[order(df$IG,df$subj.Nr, df$tp_factor),]

  # directory load scans
  write.table(dfo$scan_dir, col.names=FALSE,row.names=FALSE,quote=FALSE,
              file='scans.txt')
  # Modified SwE type - visit: tp.txt
  write.table(dfo$visit, col.names=FALSE, row.names=FALSE,quote=FALSE,
              file='tp.txt')
  # Modified SwE type - Groups: group.txt
  write.table(dfo$group, col.names=FALSE,row.names=FALSE,quote=FALSE,
              file='group.txt')
  # Subjects
  write.table(dfo$subj.ID, col.names=FALSE, row.names=FALSE,
              file='subjID.txt')
  write.table(dfo$subj.Nr, col.names=FALSE, row.names=FALSE,
              file='subjNr.txt')
  
  ## Covariates
  # nuisance covariates
  write.table(dfo$Age_BL, col.names=FALSE, row.names=FALSE, quote=FALSE,
              file='Age.txt')
  write.table(dfo$Sex, col.names=FALSE, row.names=FALSE, quote=FALSE,
              file='Sex.txt')
  write.table(dfo$logmFD, col.names=FALSE, row.names=FALSE, quote=FALSE,
              file='logmeanFD.txt')
  write.table(dfo$BMIbl, col.names=FALSE, row.names=FALSE, quote=FALSE,
              file='bmiBL.txt')
  # for covariates of interest
  write.table(dfo$IG, col.names=FALSE,row.names=FALSE, quote=FALSE,
              file='group_IG.txt')
  write.table(dfo$KG, col.names=FALSE,row.names=FALSE, quote=FALSE,
              file='group_KG.txt')
  
  # other important variables
  write.table(dfo$BMI, col.names=FALSE, row.names=FALSE, quote=FALSE,
              file='BMI.txt')
  write.table(dfo$meanFD, col.names=FALSE, row.names=FALSE, quote=FALSE,
              file='meanFD.txt')
  
  # ----------------------------------------------------------------------------
  # Different model with time as factor (instead of continous)
  # ----------------------------------------------------------------------------
  
  write.table(dfo$tp_IG, col.names=FALSE, row.names=FALSE,quote=FALSE,
            file='tp_IG.txt')
  write.table(dfo$tp_KG, col.names=FALSE, row.names=FALSE,quote=FALSE,
            file='tp_KG.txt')
  write.table(dfo$IG_fu2, col.names=FALSE,row.names=FALSE,quote=FALSE,
            file='IG_fu2.txt')
  write.table(dfo$IG_fu, col.names=FALSE,row.names=FALSE,quote=FALSE,
            file='IG_fu.txt')
  write.table(dfo$IG_bl, col.names=FALSE,row.names=FALSE,quote=FALSE,
            file='IG_bl.txt')

  write.table(dfo$KG_fu2, col.names=FALSE,row.names=FALSE,quote=FALSE,
            file='KG_fu2.txt')
  write.table(dfo$KG_fu, col.names=FALSE,row.names=FALSE,quote=FALSE,
            file='KG_fu.txt')
  write.table(dfo$KG_bl, col.names=FALSE,row.names=FALSE,quote=FALSE,
            file='KG_bl.txt')
  
  # end of computations --------------------------------------------------------
  
  # reset to original working directory
  setwd(original_wd)
}

# then save txt for all groups and return the final dataframe
get_txt_for_swe(group = "both", tp = "all", exclFD = FALSE)
get_txt_for_swe(group = "IG", tp = "all", exclFD = FALSE)
get_txt_for_swe(group = "both", tp = "BLFU", exclFD = FALSE)
get_txt_for_swe(group = "both", tp = "BL", exclFD = FALSE)

get_txt_for_swe(group = "both", tp = "all", exclFD = TRUE)
get_txt_for_swe(group = "IG", tp = "all", exclFD = TRUE)
get_txt_for_swe(group = "both", tp = "BLFU", exclFD = TRUE)
get_txt_for_swe(group = "both", tp = "BL", exclFD = TRUE)
