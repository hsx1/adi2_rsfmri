
construct_FCtables <- function() {
  
# ------------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------------

res_fc <- list()
res_fc$result_dir <- "/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis"
res_fc$BMImodel_dirs <- 
  c("noExclFD/BMI_total/brain/PCC_cc_z/bmi-age-sex-meanFD_WB-c01cl",# deact
    "noExclFD/BMI_total/brain/PCC_cc_z/bmi-age-sex_WB-c01cl")# deact
res_fc$FDBMImodel_dirs <- 
  c("noExclFD/FD_total/brain/Nacc_cc_z/bmi-fd-age-sex_WB-c02cl", # deact
    "noExclFD/FD_total/brain/Nacc_cc_z/bmi-fd-age-sex_WB-c03cl", # act
    "noExclFD/FD_total/brain/Nacc_gsr_z/bmi-fd-age-sex_WB-c03cl", # act
    "noExclFD/FD_total/brain/PCC_cc_z/bmi-fd-age-sex_WB-c01cl") # deact
res_fc$FDmodel_dirs <- 
  c("noExclFD/FD_total/brain/Nacc_cc_z/fd-age-sex_WB-c01cl", # act
    "noExclFD/FD_total/brain/Nacc_gsr_z/fd-age-sex_WB-c01cl") # act
res_fc$deact_csv <- "results_deact.csv"
res_fc$act_csv <- "results_act.csv"


# ------------------------------------------------------------------------------
# Anatomical Labelling
# ------------------------------------------------------------------------------

f <- "/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/noExclFD/result_report.txt"
result_report <- read.table(f, sep="\t",col.names = c("model","dir","k"),stringsAsFactors = FALSE)

ResultFCList <- list()
for (r in 1:nrow(result_report)) {
  f = result_report[r,"dir"]
  txt <- read.table(trimws(f), sep="\t", fill = TRUE)
  name_model <- result_report$model[r]
  
  # construct list of clusters for model
  # ClusterList <- list()
  # startentries <- which(grepl("Cluster", txt$V1))+1
  # stopentries <- which(grepl("Maximum 01", txt$V1))-1
  # for (j in 1:length(startentries)) {
  #   name_cluster <- paste('cluster',j,sep='')
  #   tmp <- txt[c(startentries[j]:stopentries[j]),c(3,5,6)]
  #   colnames(tmp) <- c("Percent", "Hem", "Label")
  #   tmp$Hem <- factor(trimws(tmp$Hem))
  #   ClusterList[[name_cluster]] <- tmp
  # }
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
    colnames(tmp) <- c("Statistics", "X","Y","Z","","Label")
    tmp[,5] <- NULL
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

# noExclFD, total sample, all tp
### first option
# f1 <- read.csv(file.path(res_fc$result_dir,res_fc$BMImodel_dirs[1],res_fc$deact_csv))
# colnr <- ncol(df1)
# df1$seed <- NULL
# df1 <- rbind(rep(NA, colnr),df1)
# df1[1:2,"seed"] <- c("\\textit{average BMI (decrease)}","PCC (cc)")
# df1$covariates <- NA
# df1$covariates[2] <- "age, sex, log mFD"
# df2 <- read.csv(file.path(res_fc$result_dir,res_fc$BMImodel_dirs[2],res_fc$deact_csv))
# df2$seed <- NULL
# df2 <- rbind(rep(NA, colnr),df2)
# df2[1:2,"seed"] <- c("\\textit{average BMI (decrease)}","PCC (cc)")
# df2$covariates <- NA
# df2$covariates[2] <- "age, sex"# BMIagesexfd_PCCcc_avgBMI_deact

### alternative (with group_row())
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
# BMIFDagesex_NACCgsr_avgFD_act 
df5 <- read.csv(file.path(res_fc$result_dir,res_fc$FDBMImodel_dirs[3],res_fc$act_csv))
df5$seed <- NULL
df5[1,"seed"] <- c("NAcc (gsr)")
df5$covariates <- NA
df5$covariates[1] <- "age, sex"
# BMIFDagesex_PCCcc_avgBMI-deact 
df6 <- read.csv(file.path(res_fc$result_dir,res_fc$FDBMImodel_dirs[4],res_fc$deact_csv))
df6$seed <- NULL
df6[1,"seed"] <- c("PCC (cc)")
df6$covariates <- NA
df6$covariates[1] <- "age, sex"

# BMIFDagesex_NACCcc_avgFD_act 
df7 <- read.csv(file.path(res_fc$result_dir,res_fc$FDmodel_dirs[1],res_fc$act_csv))
df7$seed <- NULL
df7[1,"seed"] <- c("NAcc (cc)")
df7$covariates <- NA
df7$covariates[1] <- "age, sex"
# BMIFDagesex_NACCgsr_avgFD_act
df8 <- read.csv(file.path(res_fc$result_dir,res_fc$FDmodel_dirs[2],res_fc$act_csv))
df8$seed <- NULL
df8[1,"seed"] <- c("NAcc (gsr)")
df8$covariates <- NA
df8$covariates[1] <- "age, sex"

## merge tables 
BMImodel <- rbind(df1,df2)
FDmodel <- rbind(df7,df8)
BMIFDmodel <- rbind(df3,df4,df5,df6)

ModelList <- list(BMImodel,FDmodel,BMIFDmodel)
names(ModelList) <- c("BMImodel","FDmodel","BMIFDmodel")
  
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
# Labelling
# ------------------------------------------------------------------------------

label1 <- data.table::rbindlist(
  list(
    data.frame(ResultFCList[[1]][1])[1:3, ],
    data.frame(ResultFCList[[1]][2])[1:3, ],
    data.frame(ResultFCList[[1]][3])[1:3, ],
    data.frame(ResultFCList[[2]][1])[1:3, ],
    data.frame(ResultFCList[[2]][2])[1:3, ],
    data.frame(ResultFCList[[2]][3])[1:2, ],
    data.frame(ResultFCList[[2]][4])[1:3, ]
  ),
  use.names = FALSE
)
label2 <- data.table::rbindlist(
  list(
    data.frame(ResultFCList[[7]][1])[1:3, ],
    data.frame(ResultFCList[[8]][1])[1:3, ]
    ), 
  use.names = FALSE
  )
label3 <- data.table::rbindlist(
  list(
    data.frame(ResultFCList[[3]][1])[1, ],
    data.frame(ResultFCList[[4]][1])[1:3, ],
    data.frame(ResultFCList[[5]][1])[1:3, ],
    data.frame(ResultFCList[[6]][1])[1:3, ],
    data.frame(ResultFCList[[6]][2])[1:3, ],
    data.frame(ResultFCList[[6]][3])[1:3, ]
  ),
  use.names = FALSE
)
LabelList <- list(label1,label2,label3)

for (i in 1:length(LabelList)){
  # transfer second letter of string to $hem
  hem <- stringr::str_extract(LabelList[[i]]$cluster1.Label,".{1}")
  hem <- sub('N', '', hem)
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
fn1 <- "Hem, hemisphere; L, left; R, right; MNI (Montreal Neurological Institute) coordinates of primary peak location: X = sagittal, Y = coronal, Z = axial."
fn3 <- "To identify significant clusters, we applied a cluster size threshold with p < 0.001 determined by Wild Bootstrap of 1000 samples."
fn4 <- "Three most prominent regions within a custer; more detailed description in supplement."
title_vec <- c("BMI", "FD", "BMI-FD") 

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
      caption = sprintf("Changes in functional connectivity in whole brain analysis for %s model",title_vec[i])
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
names(TableList) <- c("tab_BMImodel","tab_FDmodel","tab_BMIFDmodel")

TableList[[1]] <- TableList[[1]] %>%
  group_rows("average BMI (decrease)", 1, nrow(df1)) %>%
  group_rows("average BMI (decrease)", nrow(df1)+1, nrow(df1)+nrow(df2)) 

TableList[[2]] <- TableList[[2]] %>%
  group_rows("average log mean FD (increase)", 1, nrow(df7)) %>%
  group_rows("average log mean FD (increase)", nrow(df7)+1, nrow(df7)+nrow(df8))


TableList[[3]] <- TableList[[3]] %>%
  group_rows("change in BMI (decrease)", 1, nrow(df3)) %>%
  group_rows("average log mean FD (increase)", nrow(df3)+1, nrow(df3)+nrow(df4)) %>%
  group_rows("average log mean FD (increase)", nrow(df3)+nrow(df4)+1, nrow(df3)+nrow(df4)+nrow(df5)) %>%
  group_rows("average BMI (decrease)", nrow(df3)+nrow(df4)+nrow(df5)+1, nrow(df3)+nrow(df4)+nrow(df5)+nrow(df6))

rm(df1,df2,df3,df4,df5,df6,df7,df8)
return(TableList)
}
