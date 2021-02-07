# contact: heinrichs@cbs.mpg.de
# info: creates tables for called upon by process.R, partly as R objects, partly
# as LateX text obtained by knitr/kable. CAVE: especially FCTables requires manual
# adaptation

# Load packages -----------------------------------------------------------
library(dplyr)
library(kableExtra)


DIR_ANALYSIS <- "/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis"

# FC RESULTS -------------------------------------------------------------------
mk_FCTables <- function() {
  # previously csv files with SPM output tables need to be created, maps need
  # to be inspected with the Anatomical Labelling Toolbox and its results saved 
  # in a txt file.
  
# ------------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------------

res_fc <- list()
res_fc$result_dir <- DIR_ANALYSIS
res_fc$BMImodel_dirs <- 
  c("noExclFD/BMI_total/brain/PCC_cc_z/bmi-age-sex-meanFD_WB-c01cl",# deact
    "noExclFD/BMI_total/brain/PCC_cc_z/bmi-age-sex_WB-c01cl")# deact
res_fc$FDBMImodel_dirs <- 
  c("noExclFD/FD_total/brain/Nacc_cc_z/bmi-fd-age-sex_WB-c02cl", # deact
    "noExclFD/FD_total/brain/Nacc_cc_z/bmi-fd-age-sex_WB-c03cl", # act
    "noExclFD/FD_total/brain/Nacc_gsr_z/bmi-fd-age-sex_WB-c02cl", # deact
    "noExclFD/FD_total/brain/Nacc_gsr_z/bmi-fd-age-sex_WB-c03cl", # act
    "noExclFD/FD_total/brain/PCC_cc_z/bmi-fd-age-sex_WB-c01cl") # deact
res_fc$FDmodel_dirs <- 
  c("noExclFD/FD_total/brain/Nacc_cc_z/fd-age-sex_WB-c01cl", # act
    "noExclFD/FD_total/brain/Nacc_gsr_z/fd-age-sex_WB-c01cl") # act
res_fc$deact_csv <- "results_deact.csv"
res_fc$act_csv <- "results_act.csv"


# ------------------------------------------------------------------------------
# Anatomical Labeling
# ------------------------------------------------------------------------------

# get Anatomy tables
# manually created list of all directories of result txt - better solution needed
f <- file.path(DIR_ANALYSIS,"noExclFD/result_report.txt")
result_report <- read.table(f, sep="\t",col.names = c("model","dir","k"),stringsAsFactors = FALSE)

ResultFCList <- list()
for (r in 1:nrow(result_report)) {
  f = result_report[r,"dir"]
  txt <- read.table(trimws(f),sep="\t",fill=TRUE,col.names=paste0("V", 1:13))
  name_model <- result_report$model[r]
  
  # construct list of clusters for model
  ClusterList <- list()
  nrcluster <- which(grepl("Maximum 01", txt$V1))
  for (j in 1:length(nrcluster)) {
    name_cluster <- paste('cluster',j,sep='')
    if (j == length(nrcluster)){
      tmp <- txt[c((nrcluster[j]):nrow(txt)),c(2,4:ncol(txt))]
    }else{
      tmp <- txt[c((nrcluster[j]):nrcluster[j+1]),c(2,4:ncol(txt))]
    }
    colnames(tmp) <- c("Statistics", "X","Y","Z","","","MNIx","MNIy","MNIz","","Label")
    tmp[,c(5,6,10)] <- NULL
    tmp <- tmp[which(grepl("STAT", tmp$Statistics)),]
    ClusterList[[name_cluster]] <- tmp
  }
  # attach model with  cluster list to list of models
  ResultFCList[[name_model]] <- ClusterList
}
# data.frame(ResultFCList[[2]][3])
# ResultFCList[[2]]


# ------------------------------------------------------------------------------
# Merge tables into one for each model
# ------------------------------------------------------------------------------

# get SPM tables
# BMIagesexfd_PCCcc_avgBMI_deact 
df1 <- read.csv(file.path(res_fc$result_dir,res_fc$BMImodel_dirs[1],res_fc$deact_csv))
colnr <- ncol(df1)
df1$seed <- NULL
df1[1,"seed"] <- c("PCC (cc)")
df1$covariates <- NA
df1$covariates[1] <- "age, sex, log mFD"
# BMIagesex_PCCcc_avgBMI_deact 
df2 <- read.csv(file.path(res_fc$result_dir,res_fc$BMImodel_dirs[2],res_fc$deact_csv))
df2$seed <- NULL
df2[1,"seed"] <- c("PCC (cc)")
df2$covariates <- NA
df2$covariates[1] <- "age, sex"

# BMIFDagesex_NACCcc_cggBMI_deact 
df3 <- read.csv(file.path(res_fc$result_dir,res_fc$FDBMImodel_dirs[1],res_fc$deact_csv))
df3$seed <- NULL
df3[1,"seed"] <- c("NAcc (cc)")
df3$covariates <- NA
df3$covariates[1] <- c("age, sex")
# BMIFDagesex_NACCcc_avgFD_act
df4 <- read.csv(file.path(res_fc$result_dir,res_fc$FDBMImodel_dirs[2],res_fc$act_csv))
df4$seed <- NULL
df4[1,"seed"] <- c("NAcc (cc)")
df4$covariates <- NA
df4$covariates[1] <- "age, sex"
# BMIFDagesex_NACCgsr_cggBMI_deact 
df5 <- read.csv(file.path(res_fc$result_dir,res_fc$FDBMImodel_dirs[3],res_fc$deact_csv))
df5$seed <- NULL
df5[1,"seed"] <- c("NAcc (gsr)")
df5$covariates <- NA
df5$covariates[1] <- "age, sex"
# BMIFDagesex_NACCgsr_avgFD_act 
df6 <- read.csv(file.path(res_fc$result_dir,res_fc$FDBMImodel_dirs[4],res_fc$act_csv))
df6$seed <- NULL
df6[1,"seed"] <- c("NAcc (gsr)")
df6$covariates <- NA
df6$covariates[1] <- "age, sex"
# BMIFDagesex_PCCcc_avgBMI-deact 
df7 <- read.csv(file.path(res_fc$result_dir,res_fc$FDBMImodel_dirs[5],res_fc$deact_csv))
df7$seed <- NULL
df7[1,"seed"] <- c("PCC (cc)")
df7$covariates <- NA
df7$covariates[1] <- "age, sex"

# BMIFDagesex_NACCcc_avgFD_act 
df8 <- read.csv(file.path(res_fc$result_dir,res_fc$FDmodel_dirs[1],res_fc$act_csv))
df8$seed <- NULL
df8[1,"seed"] <- c("NAcc (cc)")
df8$covariates <- NA
df8$covariates[1] <- "age, sex"
# BMIFDagesex_NACCgsr_avgFD_act
df9 <- read.csv(file.path(res_fc$result_dir,res_fc$FDmodel_dirs[2],res_fc$act_csv))
df9$seed <- NULL
df9[1,"seed"] <- c("NAcc (gsr)")
df9$covariates <- NA
df9$covariates[1] <- "age, sex"

## merge tables 
BMImodel <- rbind(df2,df1)
BMIFDmodel <- rbind(df3,df4,df5,df6,df7)
FDmodel <- rbind(df8,df9)

ModelList <- list(BMImodel,BMIFDmodel,FDmodel)
names(ModelList) <- c("BMImodel","BMIFDmodel","FDmodel")
  
for (i in 1:length(ModelList)) {
  ModelList[[i]]$anat <- NA
  ModelList[[i]]$hem <- NA
  # select columns
  ModelList[[i]] <-
    ModelList[[i]][, c(
      "seed",
      "covariates",
      "cluster_p.FWE.corr.",
      "cluster_equivk",
      "cluster_equivkZ",
      "peak_p.FWE.corr.",
      "X_x",
      "y",
      "z..mm.",
      "hem",
      "anat"
    )]
  # name columns
  colnames(ModelList[[i]]) <-
    c(
      "Seed",
      "Covariates",
      "FWE-corr. P",
      "cluster size",
      "Z score", # of custer size
      "FWE-corr. P",
      "X",
      "Y",
      "Z",
      "Hem",
      "Anatomical region²"
    )
}

# only keep cluster information
# for (i in 1:length(ModelList)){
#   # only keep cluster statistics and peak local maximum
#   ModelList[[i]] <- ModelList[[i]][!is.na(ModelList[[i]]$Seed) | !is.na(ModelList[[i]]$`clusterwise FWE-corr. P`),]
# }

# ------------------------------------------------------------------------------
# Labeling
# ------------------------------------------------------------------------------

# wird alles überarbeitet
# BMI
label1 <- data.table::rbindlist(
  list(
    # for cluster of df1, df2
    # avgBMI - age,sex
    data.frame(ResultFCList[[2]][1])[1:3, ] %>% select(last_col()), # cl1
    data.frame(ResultFCList[[2]][2])[1:2, ] %>% select(last_col()),
    data.frame(ResultFCList[[2]][3])[1:2, ] %>% select(last_col()),
    data.frame(ResultFCList[[2]][4])[1:2, ] %>% select(last_col()),
    # avgBMI - age,sex, mFD
    data.frame(ResultFCList[[1]][1])[1:3, ] %>% select(last_col()), # cl1
    data.frame(ResultFCList[[1]][2])[1:2, ] %>% select(last_col()), # cl2
    data.frame(ResultFCList[[1]][3])[1:2, ] %>% select(last_col()), # ...
    data.frame(ResultFCList[[1]][4])[1:2, ] %>% select(last_col())
  ),
  use.names = FALSE
)
# BMI FD
label2 <- data.table::rbindlist(
  list(
    # cgnBMI - NAcc cc age,sex
    data.frame(ResultFCList[[3]][1])[1, ] %>% select(last_col()),
    # avgFD - NAcc cc age,sex
    data.frame(ResultFCList[[4]][1])[1:3, ] %>% select(last_col()),
    # cgnBMI - NAcc gsr age,sex
    data.frame(ResultFCList[[5]][1])[1:2, ] %>% select(last_col()),
    # avgFD - NAcc gsr age,sex
    data.frame(ResultFCList[[6]][1])[1:3, ] %>% select(last_col()),
    # avgBMI - PCC cc age,sex
    data.frame(ResultFCList[[7]][1])[1:3, ] %>% select(last_col()),
    data.frame(ResultFCList[[7]][2])[1:2, ] %>% select(last_col()),
    data.frame(ResultFCList[[7]][3])[1:2, ] %>% select(last_col()),
    data.frame(ResultFCList[[7]][4])[1:2, ] %>% select(last_col())
  ),
  use.names = FALSE
)
# FD
label3 <- data.table::rbindlist(
  list(
    data.frame(ResultFCList[[8]][1])[1:3, ] %>% select(last_col()),
    data.frame(ResultFCList[[9]][1])[1:3, ] %>% select(last_col())
  ), 
  use.names = FALSE
)
LabelList <- list(label1,label2,label3)

for (i in 1:length(LabelList)){
  # transfer second letter of string to $hem
  hem <- stringr::str_extract(LabelList[[i]]$cluster1.Label,".{1}")
  hem <- sub('[^RL]', '', hem)
  hem[is.na(hem)] <- ''
  LabelList[[i]]$Hem <- hem
  # remove fist three letter from all strings
  LabelList[[i]]$cluster1.Label <- sub('..', '', LabelList[[i]]$cluster1.Label)
  LabelList[[i]]$cluster1.Label <- sub('A', '-', LabelList[[i]]$cluster1.Label)
}

ModelList[[1]]$`Anatomical region²` <- LabelList[[1]]$cluster1.Label
ModelList[[1]]$Hem <- LabelList[[1]]$Hem
ModelList[[2]]$`Anatomical region²` <- LabelList[[2]]$cluster1.Label
ModelList[[2]]$Hem <- LabelList[[2]]$Hem
ModelList[[3]]$`Anatomical region²` <- LabelList[[3]]$cluster1.Label
ModelList[[3]]$Hem <- LabelList[[3]]$Hem


# ------------------------------------------------------------------------------
# Construct Latex Tables
# ------------------------------------------------------------------------------

# hide NA in tables
fn1 <- "Hem, hemisphere; L, left; R, right; FWE-corr., family-wise error corrected; MNI (Montreal Neurological Institute) coordinates of primary peak location: X sagittal; Y, coronal; Z, axial; cc, preprocessing with AROMA-ICA + CompCor; gsr, preprocessing with AROMA-ICA + CompCor + GSR."
fn3 <- "To identify significant clusters, we applied a cluster size threshold with p < 0.001 determined by Wild Bootstrap of 1000 samples."
fn4 <- "Connectivity with maximum three voxels that mark local maxima within the respective custer; more detailed description of anatomical regions that are assigned to overall clusters and and corresponding probability in Supplementary."
title_vec <- c("BMI", "BMI-FD", "FD") 

tabnames <- c("tableBMImodel","tableBMIFDmodel","tableFDmodel")
TableList <- list()
for (i in 1:length(ModelList)){
  TableList[[i]] <-
    ModelList[[i]] %>%
    knitr::kable(
      escape = FALSE,     # use font spec of latex with kableExtra
      col.names = colnames(ModelList[[i]]),
      row.names = FALSE,
      format = "latex",
      booktabs = T,
      linesep = "",      # disable "\\addlinespace" at every 5th line
      label = tabnames[i],
      caption = sprintf("Changes in functional connectivity in whole brain analysis for %s model",title_vec[i]),
    ) %>%
    kable_styling(latex_options = c('scale_down')) %>% # "scale_down"
    #column_spec(c(11), width = "3.2cm") %>%
    #column_spec(c(2:6), width = "2.5cm") %>%
    #column_spec(c(4,9,10,11), width = "0.5cm") %>%
    add_header_above(header = c(" " = 6,"MNI coordinates" = 3)) %>%
    add_header_above(header = c(" "," ","clusterwise¹" = 3,"voxelwise  at local maximum" = 6)) %>%
    add_footnote(fn1, notation = "none") %>%
    add_footnote(c(fn3, fn4), notation = "number") %>%
    kableExtra::landscape()
}
names(TableList) <- tabnames

TableList[[1]] <- TableList[[1]] %>%
  kableExtra::group_rows("average BMI (decrease)", 1, nrow(df2)) %>%
  kableExtra::group_rows("average BMI (decrease)", nrow(df2)+1, nrow(df2)+nrow(df1)) 

TableList[[2]] <- TableList[[2]] %>%
  kableExtra::group_rows("change in BMI (decrease)", 1, nrow(df3)) %>%
  kableExtra::group_rows("average log mean FD (increase)", nrow(df3)+1, nrow(df3)+nrow(df4)) %>%
  kableExtra::group_rows("change in BMI (decrease)", nrow(df3)+nrow(df4)+1, nrow(df3)+nrow(df4)+nrow(df5)) %>%
  kableExtra::group_rows("average log mean FD (increase)", nrow(df3)+nrow(df4)+nrow(df5)+1, nrow(df3)+nrow(df4)+nrow(df5)+nrow(df6)) %>%
  kableExtra::group_rows("average BMI (decrease)", nrow(df3)+nrow(df4)+nrow(df5)+nrow(df6)+1, nrow(df3)+nrow(df4)+nrow(df5)+nrow(df6)+nrow(df7))

TableList[[3]] <- TableList[[3]] %>%
  kableExtra::group_rows("average log mean FD (increase)", 1, nrow(df8)) %>%
  kableExtra::group_rows("average log mean FD (increase)", nrow(df8)+1, nrow(df8)+nrow(df9))

rm(df1,df2,df3,df4,df5,df6,df7,df8)
return(TableList)
}

# SAMPLE DESCRIPTION -----------------------------------------------------------

mk_SampleTable <- function(final) {
  # completed time points
  freqtable <- as.data.frame(table(final$subj.ID))
  tpfreq <- table(freqtable$Freq)
  subj_df <- final[match(unique(final$subj.ID), final$subj.ID),]
  
  # overview distribution of data points
  tp1 <- as.character(freqtable$Var1[freqtable$Freq == 1])
  tp2 <- as.character(freqtable$Var1[freqtable$Freq == 2])
  tp3 <- as.character(freqtable$Var1[freqtable$Freq == 3])
  tp1_tab <- table(subset(final, subj.ID %in% tp1)$tp, subset(final, subj.ID %in% tp1)$condition)
  tp2_tab <- table(subset(final, subj.ID %in% tp2)$tp, subset(final, subj.ID %in% tp2)$condition)[c(1,3),]
  tp3_tab <- table(subset(final, subj.ID %in% tp3)$tp, subset(final, subj.ID %in% tp3)$condition)[1,]
  total_datapoints <- c(sum(final$group==1),sum(final$group==2))
  total_subjects <- plyr::count(subj_df$group)$freq
  tab_sample <- data.frame(rbind(tp1_tab, tp2_tab, tp3_tab, total_subjects, total_datapoints))
  rm(tp1,tp2,tp3,tp1_tab, tp2_tab, tp3_tab, total_subjects, total_datapoints)
  rownames(tab_sample) <- c("count: only 0","count: only 6","count: only 12","count: 0 and 6","count: 6 and 12","count: complete data","total number of subjects","total data points")
  colnames(tab_sample) <- c("BARS","NBARS")
  return(tab_sample)
}




# Detailed FC Labels (supplement) -----------------------------------------

mk_DetailedLabelTab <- function(){
  # Contains two List, one with all Tables $EntryList, and one with all table 
  # captions $CaptionList for alle Tables in the same order.
  DIR_ANALYSIS <- "/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis"
  df <- file.path(DIR_ANALYSIS,"noExclFD/result_report.txt")
  AnatomyResults <- read.table(df, sep="\t",col.names = c("label","dir","k"),stringsAsFactors = FALSE)
  
  
  EntryList <- list()
  
  for (r in 1:nrow(AnatomyResults)) {
    
    file = AnatomyResults[r,"dir"]
    txt_cols <- c("VoxelCount","equals","PercClusterVolumeAssignedTo","in","Hem","Region","PercentOfArea","V8","V9","V10","V11","V12","V13")
    txt <- read.table(trimws(file),sep="\t",fill=TRUE,col.names=txt_cols) # paste0("V", 1:13)
    
    
    
    txt <- txt[sort(c(which(grepl("Cluster", txt$VoxelCount)),which(grepl("voxel", txt$equals)))),]
    txt <- txt[,c("VoxelCount","PercClusterVolumeAssignedTo", "Hem","Region","PercentOfArea")]
    
    idxClusterStart <- which(grepl("Cluster", txt$VoxelCount))
    
    #txt$VoxelCount <- as.character(txt$VoxelCount)
    txt <- cbind(descr = NA, txt)
    txt$descr[idxClusterStart] <- stringr::str_remove(string = as.character(txt$VoxelCount[idxClusterStart]), pattern = "\\:[^:]*$")
    txt$VoxelCount[idxClusterStart] <- NA
    
    # attach model with  cluster list to list of models
    name_model <- AnatomyResults$label[r]
    EntryList[[name_model]] <- txt
  }
  
  
  # Construct table title ---------------------------------------------------
  
  AnatomyResults <-
    tidyr::separate(
      data = AnatomyResults,
      col = "label",
      into = c("label", "roi", "variable", "association"),
      sep = "_"
    )
  
  AnatomyResults$start <-
    "Anatomical labels for clusters functionally connected to the "
  AnatomyResults <-
    tidyr::extract(
      data = AnatomyResults,
      col = "roi",
      into = "denoising",
      regex = "([a-z]+)",
      remove = FALSE
    ) %>%
    tidyr::extract(col = "roi", into = "roi", regex = "([A-Z]+)")
  
  AnatomyResults$denoising <- AnatomyResults$denoising %>%
    stringr::str_replace_all('cc', 'AROMA+CC') %>%
    stringr::str_replace_all('gsr', 'AROMA+CC+GSR')
  
  AnatomyResults$association <- AnatomyResults$association %>%
    stringr::str_replace_all('deact', 'negative') %>%
    stringr::str_replace_all('act', 'positive')
  
  AnatomyResults$variable <- AnatomyResults$variable %>%
    stringr::str_replace_all('cng', 'change in ') %>%
    stringr::str_replace_all('avg', 'average ') %>%
    stringr::str_replace_all('FD', 'mFD')
  
  AnatomyResults <-
    tidyr::extract(
      data = AnatomyResults,
      col = "label",
      into = "adjust",
      regex = "([a-z]+)",
      remove = FALSE
    ) %>%
    tidyr::extract(col = "label", into = "model", regex = "([A-Z]+)")
  
  AnatomyResults$model <- AnatomyResults$model %>%
    stringr::str_replace_all('BMIFD', 'BMI-FD')
  
  AnatomyResults$adjust <- AnatomyResults$adjust %>%
    stringr::str_replace_all('agesexfd', 'age, sex and mFD') %>%
    stringr::str_replace_all('agesex', 'age and sex')
  
  AnatomyResults$roi <- AnatomyResults$roi %>%
    stringr::str_replace_all('NACC', 'NAcc')
  
  CaptionList <- with(
    AnatomyResults,
    paste(
      "Anatomical labels for clusters functionally connected to the",
      roi, 
      "showing a",
      association,
      "association with",
      variable,
      "adjusted for",
      adjust,
      "in the",
      model,
      "model on",
      denoising,
      "preprocessed data."
    )
  )
  
  # knitr tables ------------------------------------------------------------
  
  txt_colname <- kableExtra::linebreak(c("Cluster","Number of voxels\n in cluster","\\% of cluster volume\nassigned", "Hemisphere","Area","\\% of area overlap\nwith cluster"))
  
  ToolboxOutputList <- list(EntryList, CaptionList)
  names(ToolboxOutputList) <- c("EntryList", "CaptionList")
  return(ToolboxOutputList)
}
