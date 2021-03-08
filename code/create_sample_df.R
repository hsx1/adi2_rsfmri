library(dplyr) # version 1.0.2
library(tidyr) # version 1.1.2

CONNMAPS_PATH <- "/data/pt_02161/Results/Project2_resting_state/connectivity/calc_DMN_reward_seed_connectivity"
DF_PATH <- "/data/p_02161/ADI_studie/metadata/final_sample_MRI_QA_info.csv"
ADDITIONALINFO_FILE <-"/data/pt_02161/Data/metadata/Prehn_metadata/ADI_Rekrutierung_noPW.xls"

create_sample_df <- function(group = "both", tp = "all", exclFD = FALSE, out = "final"){
  # group = IG/KG/both
  # tp = BL/FU/FU2/BLFU/FUFU2/all
  
  # ----------------------------------------------------------------------------
  # primary definitions
  # ----------------------------------------------------------------------------
  
  # import absolute path
  abs_path = read.csv("../abs_path.csv",
                      header = FALSE,
                      stringsAsFactors = FALSE)
  
  # ----------------------------------------------------------------------------
  # functions
  # ----------------------------------------------------------------------------
  
  
  apply_excl_criteria <-
    function(full = full,
             mri = mri_files,
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
          "comorbidities",
          "interv",
          "intervention_de",
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
  
  mri_files = data.frame(list.files(CONNMAPS_PATH,
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
  
  full = read.csv(DF_PATH) # CAVE: info file elsewhere !!!
  full = full[, c(
    "subj.ID",
    "condition",
    "tp",
    "Age_BL",
    "Sex",
    "meanFD",
    "maxFD",
    "BMI",
    "Final_Score",
    "Excessive_motion",
    "Exclude",
    "Check_comments"
  )]
  
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
    seed = 8745
  )
  imp$predictorMatrix
  summary(imp)
  
  # compute mean from imputed values in complete datasets 
  com <- suppressMessages(complete(imp, "repeated"))
  com <- cbind(tmp,com[,(ncol(com)-(dfs-1)):ncol(com)])
  meanval <- apply(com[,(ncol(tmp)+1):ncol(com)],1,mean)
  sdval <- apply(com[,(ncol(tmp)+1):ncol(com)],1,sd)
  rm(com, ini, imp, tmp, pred)
  
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
  
  
  # add surgery & comorbidities ------------------------------------------------
  
  addit_info <- data.frame(t(readxl::read_excel(ADDITIONALINFO_FILE,
                                                sheet = 2,
                                                col_names = FALSE)
                             ), 
                           stringsAsFactors = FALSE
                           )
  colnames(addit_info) <- addit_info[1,]
  addit_info <- addit_info[2:nrow(addit_info),c(1:38)]
  addit_info$subj.ID <- gsub("-", "", addit_info$Code)
  comorbidities <- addit_info %>%
    dplyr::select(c("Adipositas durch Kalorienzufuhr, [BMI] 40 und mehr":"Migräne")) %>%
    sapply(as.numeric, na.rm=T) %>%
    data.frame() %>%
    sapply(as.logical, na.rm=T) %>%
    data.frame() %>%
    rowSums(., na.rm = T) 
  addit_info$comorbidities <- comorbidities
  rm(comorbidities)
  
  addit_bl <- addit_info
  addit_bl$tp <- "bl"
  addit_fu <- addit_info
  addit_fu$tp <- "fu"
  addit_fu2 <- addit_info
  addit_fu2$tp <- "fu2"
  addit_info <- rbind(addit_bl,addit_fu,addit_fu2)
  rm(addit_bl,addit_fu,addit_fu2)
  addit_info$intervention_de <- addit_info$Operationsverfahren
  addit_info$tp <- as.factor(addit_info$tp)
  addit_info$subj.ID <- as.factor(addit_info$subj.ID)

  
  full_surg_comorbid <- merge(full, addit_info, by=c("subj.ID","tp"), all = TRUE)
  rm(addit_info)
  # check if KG have 
  t(table(full_surg_comorbid$condition, full_surg_comorbid$intervention_de))
  full_surg_comorbid$interv <- NA
  full_surg_comorbid$interv[stringr::str_detect(full_surg_comorbid$intervention_de, stringr::regex("resektion|schlauch",ignore_case = T)) & !is.na(full_surg_comorbid$intervention_de)] <- "VSG"
  full_surg_comorbid$interv[stringr::str_detect(full_surg_comorbid$intervention_de, stringr::regex("band",ignore_case = T)) & !is.na(full_surg_comorbid$intervention_de)] <- "GB"
  full_surg_comorbid$interv[stringr::str_detect(full_surg_comorbid$intervention_de, stringr::regex("bypass",ignore_case = T)) & !is.na(full_surg_comorbid$intervention_de)] <- "RYGB"
  full_surg_comorbid$interv[stringr::str_detect(full_surg_comorbid$intervention_de, "KG") & !is.na(full_surg_comorbid$intervention_de)] <- "no"
  full_surg_comorbid[,c("interv", "intervention_de")]
  full_surg_comorbid$intervention 
  
  
  # Exclusion of participants --------------------------------------------------
  final_sample <-
    apply_excl_criteria(full_surg_comorbid, mri_files, exclFD)
  
  df = suppressMessages(plyr::join(mri_files, final_sample))
  df = df[!is.na(df$include), ]
  rm(apply_excl_criteria, mri_files, final_sample, i)
  
  
  # preparation of variables ---------------------------------------------------
  df$logmFD <- log10(df$meanFD)
  
  df$group[df$condition == "IG"] = 1
  df$group[df$condition == "KG"] = 2
  df$group_factor = as.factor(df$group)
  levels(df$group_factor) = c("IG", "KG")
  
  df$tp = as.factor(df$tp)
  levels(df$tp) = c("bl", "fu", "fu2")
  
  df$visit = df$tp
  levels(df$visit) = c(1, 2, 3)
  
  df$Sex = as.factor(df$Sex)
  levels(df$Sex) = c(-1, 1)
  
  
  # Modify Design matrix for subsamples ----------------------------------------
  
  # Model specification incl. design matrix for different options
  # selection of groups
  if (group == "both") {
  # total: both groups = original dataframe
  } else if (group == "IG") {
    df <- df[df$condition == "IG", ]
  } else if (group == "KG") {
    df <- df[df$condition == "KG", ]
  }
  
  # selection of time points
  if (tp == "all") {
    df <- df
  } else if (tp == "BL") {
    df <- df[df$tp == "bl", ]
  } else if (tp == "FU") {
    df <- df[df$tp == "fu", ]
  } else if (tp == "FU2") {
    df <- df[df$tp == "fu2", ]
  } else if (tp == "BLFU") {
    df <- df[df$tp != "fu2", ]
  } else if (tp == "FUFU2") {
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
  
  
  # end of computations --------------------------------------------------------
  if (out == "final") {
    return(df)
  } else if (out == "full") {
    return(full_surg_comorbid)
  }
}
