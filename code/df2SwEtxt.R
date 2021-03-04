library(dplyr) # version 1.0.2
library(tidyr) # version 1.1.2


ROOT = "/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /code_and_manuscript/"
ConnectMapsPath <-
  "/data/pt_02161/Results/Project2_resting_state/connectivity/calc_DMN_reward_seed_connectivity"
dfPath <-
  "../../../../../../p_02161/ADI_studie/metadata/final_sample_MRI_QA_info.csv"
setwd(ROOT)

get_txt_for_swe <- function(group = "both", tp = "all", exclFD = FALSE){
  # group = IG/KG/both
  # tp = BL/FU/FU2/BLFU/FUFU2/all
  
  # ----------------------------------------------------------------------------
  # primary definitions
  # ----------------------------------------------------------------------------
  
  original_wd = getwd()
  # import absolute path
  abs_path = read.csv("abs_path.csv",
                      header = FALSE,
                      stringsAsFactors = FALSE)
  # set working directory
  parentdir = file.path(getwd(), "../SwE_files", fsep = .Platform$file.sep)[1]
  setwd(parentdir)
  parentdir = getwd()
  
  
  # ----------------------------------------------------------------------------
  # second level functions
  # ----------------------------------------------------------------------------
  
  apply_excl_criteria <-
    function(full = full,
             mri = mri_sample,
             exclFD = FALSE) {
      # exclusion due to problems in preprocessing
      cat("N =", nrow(full), "data points.\n")
      full$exclude_prep <- full$Exclude
      full$exclude_prep[is.na(full$exclude_prep)] = FALSE
      mri_sample <-
        merge(mri_files, full[!is.na(full$condition), c(
          "subj.ID",
          "condition",
          "tp",
          "Age_BL",
          "Sex",
          "meanFD",
          "BMI",
          "BMI_BL",
          "BMI_BLi",
          "exclude_prep"
        )], by = c("subj.ID", "tp"))
      cat("n =", nrow(mri_sample), "fMRI data points, ")
      
      mri_sample[mri_sample$exclude_prep == TRUE, ] # ADI063_bl ADI063_fu ADI116_bl ADI116_fu ADI116_fu2
      cat(
        "of which",
        nrow(mri_sample[mri_sample$exclude_prep == TRUE, ]),
        "data points were excluded due to problems in preprocessing.\n"
      )
      sample_after_1stExcl <-
        mri_sample[mri_sample$exclude_prep == FALSE, ]
      sample_after_1stExcl$exclude_prep <- NULL
      
      if (exclFD == TRUE) {
        # exclusion of 10% of the worst mean FD
        thres10 = sort(sample_after_1stExcl$meanFD, decreasing = TRUE)[round(nrow(sample_after_1stExcl) /
                                                                               10)]
        sample_after_1stExcl$exclude_fd = (sample_after_1stExcl$meanFD >= thres10)
        cat(
          "Note that",
          nrow(sample_after_1stExcl[sample_after_1stExcl$exclude_fd == TRUE,]),
          "data points were excluded due to extensive head motion.\n"
        )
        
        final_sample = sample_after_1stExcl[sample_after_1stExcl$exclude_fd ==
                                              FALSE, ]
        sample_after_1stExcl$exclude_fd <- NULL
      } else if (exclFD == FALSE) {
        final_sample = sample_after_1stExcl
      }
      cat(
        "The remaining sample comprises",
        nrow(final_sample),
        "data points across time for all groups.\n"
      )
      # all remaining data points should be included
      final_sample$include <- TRUE
      return(final_sample)
    }
  
  # ----------------------------------------------------------------------------
  # script
  # ----------------------------------------------------------------------------
  
  # participants with mri data (txt file with path to .nii file) ---------------
  
  mri_files = data.frame(list.files(ConnectMapsPath,
                                    full.names = TRUE))
  # split info from path in txt
  for (i in 1:nrow(mri_files)) {
    tmp = strsplit(toString(mri_files[i, 1]), '/')[[1]][8]
    mri_files[i, "scan_dir"] = paste0(strsplit(toString(mri_files[i, 1]), '/')[[1]][1:8],
                                      sep = "/",
                                      collapse = "")
    mri_files[i, "tp"] = strsplit(toString(tmp), '_')[[1]][2]
    mri_files[i, "subj.ID"] = strsplit(toString(tmp), '_')[[1]][1]
    mri_files[i, "subj.Nr"] = as.numeric(substr(mri_files$subj.ID[i], 4, 6))
  }
  mri_files[, 1] <- NULL
  # exclude other time points than bl fu or fu2
  mri_files <-
    mri_files[stringr::str_detect(mri_files$tp, "^(bl|fu|fu2)$"), ]
  
  # load info from full sample -------------------------------------------------
  
  full = read.csv(dfPath) # CAVE: info file elsewhere !!!
  full = full[, c("subj.ID",
                  "condition",
                  "tp",
                  "Age_BL",
                  "Sex",
                  "meanFD",
                  "BMI",
                  "Exclude")]
  
  # add variable BMI_BL for all time points
  full$BMI_BL <- NA
  full$BMI_BL[full$tp == "bl"] <-
    full$BMI[full$tp == "bl"]
  for (i in 1:nrow(full)) {
    # check if bl exist for subject
    if (any(full$tp[full$subj.ID == full$subj.ID[i]] == "bl")) {
      # replace with value of same subject at bl
      full$BMI_BL[i] <-
        full$BMI_BL[full$subj.ID == full$subj.ID[i] &
                      full$tp == "bl"]
    }
  }
  
  # missingness -------------------------------------------------------------
  
  VIM::barMiss(full[ ,c("condition","BMI")])
  VIM::histMiss(full[,c("BMI_BL","BMI")])
  
  
  # mice for imputation in multilevel data ---------------------------------------
  
  # see: https://stackoverflow.com/questions/47950304/random-effects-in-longitudinal-multilevel-imputation-models-using-mice
  # imputation based on following variables
  tmp <- full[,c("subj.ID","condition", "tp", "Age_BL", "BMI")]
  tmp$subj.ID <- as.integer(tmp$subj.ID)
  
  # set imputation method
  impmethod <- c("", "","", "", "2l.pan") # why 2l.pan?
  names(impmethod) <- colnames(tmp)
  
  # set default predictor matrix
  ini <- mice::mice(tmp, maxit=0)
  
  # create imputation model for BMI, predictor matrix (following 'type')
  pred <- ini$predictorMatrix
  pred["BMI",] <- c(-2, 1, 2, 1, 0) 
  # (cluster variable, fixef with randint, fixef with randint + randslo, dependent)
  
  dfs = 50
  iter = 10
  imp <- mice::mice(
    tmp,
    method = impmethod,
    pred = pred,
    maxit = iter,
    m = dfs,
    seed = 8745 # DO NOT CHANGE (same in create_sample_df.R)
  )
  imp$predictorMatrix
  summary(imp)
  
  # compute mean from imputed values in complete datasets 
  com <- complete(imp, "repeated")
  com <- cbind(tmp,com[,(ncol(com)-(dfs-1)):ncol(com)])
  meanval <- apply(com[,(ncol(tmp)+1):ncol(com)],1,mean)
  sdval <- apply(com[,(ncol(tmp)+1):ncol(com)],1,sd)
  
  # impute data only for baseline (bl) and construct BMI_BLi
  full$BMI_BLi <- meanval # imputed data for all tp
  full$BMI_BLi[full$tp != "bl"] <- NA
  full$BMI_BLi[full$tp == "bl"] <-
    meanval[full$tp == "bl"]
  for (i in 1:nrow(full)) {
    # check if bl exist for subject
    if (any(full$tp[full$subj.ID == full$subj.ID[i]] == "bl")) {
      # replace with value of same subject at bl
      full$BMI_BLi[i] <-
        full$BMI_BLi[full$subj.ID == full$subj.ID[i] &
                       full$tp == "bl"]
    }
  }
  full[,c("BMI_BL","BMI_BLi")]
  
  # Exclusion of participants without mri data ---------------------------------
  final_sample <-
    apply_excl_criteria(full, mri_files, exclFD)
  
  df = suppressMessages(plyr::join(mri_files, final_sample))
  df = df[!is.na(df$include), ]
  rm(apply_excl_criteria,
     mri_files,
     full,
     final_sample,
     tmp,
     i)
  
  
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
  
  df$tp = as.factor(df$tp)
  df$visit = df$tp
  levels(df$visit) = c(1, 2, 3)
  df$tp_cov = df$tp
  levels(df$tp_cov) = c(-1, 0, 1)
  
  df$Sex = as.factor(df$Sex)
  levels(df$Sex) = c(-1, 1)
  
  # linear modeling of time for each group
  df$tp_IG <- df$tp_cov
  df[df$KG == 1, "tp_IG"] = 0
  df$tp_KG <- df$tp_cov
  df[df$IG == 1, "tp_KG"] = 0
  
  # model with time as factor
  df[df$tp == "fu2" & df$IG == 1, "IG_fu2"] = 1
  df[df$tp == "fu" & df$IG == 1, "IG_fu"] = 1
  df[df$tp == "bl" & df$IG == 1, "IG_bl"] = 1
  df[is.na(df$IG_fu2), "IG_fu2"] = 0
  df[is.na(df$IG_fu), "IG_fu"] = 0
  df[is.na(df$IG_bl), "IG_bl"] = 0
  
  df[df$tp == "fu2" & df$KG == 1, "KG_fu2"] = 1
  df[df$tp == "fu" & df$KG == 1, "KG_fu"] = 1
  df[df$tp == "bl" & df$KG == 1, "KG_bl"] = 1
  df[is.na(df$KG_fu2), "KG_fu2"] = 0
  df[is.na(df$KG_fu), "KG_fu"] = 0
  df[is.na(df$KG_bl), "KG_bl"] = 0
  
  levels(df$tp) = c("bl", "fu", "fu2")
  
  # Modify Design matrix for subsamples -------------------------------------
  
  # Model specification incl. design matrix for different options
  output_dir <- parentdir
  if (exclFD == FALSE) {
    output_dir <- file.path(parentdir, "noExclFD")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
  } else if (exclFD == TRUE) {
    output_dir <- file.path(parentdir, "ExclFD")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
  }
  
  # selection of groups
  if (group == "both") {
    output_dir <- file.path(output_dir, "both")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
  } else if (group == "IG") {
    output_dir <- file.path(output_dir, "IG")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    df <- df[df$condition == "IG", ]
  } else if (group == "KG") {
    output_dir <- file.path(output_dir, "KG")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    df <- df[df$condition == "KG", ]
  }
  
  # selection of time points
  
  if (tp == "all") {
    output_dir <- file.path(output_dir, "total")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    setwd(output_dir)
    df <- df
  } else if (tp == "BL") {
    output_dir <- file.path(output_dir, "BL")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    setwd(output_dir)
    df <- df[df$tp == "bl", ]
  } else if (tp == "FU") {
    output_dir <- file.path(output_dir, "FU")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    setwd(output_dir)
    df <- df[df$tp == "fu", ]
  } else if (tp == "FU2") {
    output_dir <- file.path(output_dir, "FU2")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    setwd(output_dir)
    df <- df[df$tp == "fu2", ]
  } else if (tp == "BLFU") {
    output_dir <- file.path(output_dir, "BLFU")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    setwd(output_dir)
    df <- df[df$tp != "fu2", ]
  } else if (tp == "FUFU2") {
    output_dir <- file.path(output_dir, "FUFU2")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    setwd(output_dir)
    df <- df[df$tp != "fu2", ]
  }
  cat(
    "Data from",
    group,
    "group(s) and",
    tp,
    "measurements were selected to be kept in the data frame, therefore n =",
    nrow(df),
    ".\n"
  )
  cat("Output directory of txt files:", output_dir)
  
  
  # order dataframe by group / subject / time point
  df$tp_factor <- as.numeric(df$tp)
  dfo <- df[order(df$IG, df$subj.Nr, df$tp_factor), ]
  
  # directory load scans
  write.table(
    dfo$scan_dir,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'scans.txt'
  )
  # Modified SwE type - visit: tp.txt
  write.table(
    dfo$visit,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'tp.txt'
  )
  # Modified SwE type - Groups: group.txt
  write.table(
    dfo$group,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'group.txt'
  )
  # Subjects
  write.table(
    dfo$subj.ID,
    col.names = FALSE,
    row.names = FALSE,
    file = 'subjID.txt'
  )
  write.table(
    dfo$subj.Nr,
    col.names = FALSE,
    row.names = FALSE,
    file = 'subjNr.txt'
  )
  
  ## Covariates
  # nuisance covariates
  write.table(
    dfo$Age_BL,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'Age.txt'
  )
  write.table(
    dfo$Sex,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'Sex.txt'
  )
  write.table(
    dfo$logmFD,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'logmeanFD.txt'
  )
  # imputed baseline BMI
  write.table(
    dfo$BMI_BLi,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'BMI_BLi.txt'
  )
  # for covariates of interest
  write.table(
    dfo$IG,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'group_IG.txt'
  )
  write.table(
    dfo$KG,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'group_KG.txt'
  )
  
  # other important variables
  write.table(
    dfo$BMI,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'BMI.txt'
  )
  write.table(
    dfo$meanFD,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'meanFD.txt'
  )
  
  # ----------------------------------------------------------------------------
  # Different model with time as factor (instead of continuous)
  # ----------------------------------------------------------------------------
  
  write.table(
    dfo$tp_IG,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'tp_IG.txt'
  )
  write.table(
    dfo$tp_KG,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'tp_KG.txt'
  )
  write.table(
    dfo$IG_fu2,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'IG_fu2.txt'
  )
  write.table(
    dfo$IG_fu,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'IG_fu.txt'
  )
  write.table(
    dfo$IG_bl,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'IG_bl.txt'
  )
  
  write.table(
    dfo$KG_fu2,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'KG_fu2.txt'
  )
  write.table(
    dfo$KG_fu,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'KG_fu.txt'
  )
  write.table(
    dfo$KG_bl,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    file = 'KG_bl.txt'
  )
  
  # end of computations --------------------------------------------------------
  
  # reset to original working directory
  setwd(original_wd)
  return(dfo)
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
