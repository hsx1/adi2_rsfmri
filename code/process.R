# processes data to report later in manuscript .Rmd file


# Get functions -----------------------------------------------------------
setwd("/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /code_and_manuscript/code/")
source("figures.R")
source("tables.R")
source("create_sample_df.R")


# Load data ---------------------------------------------------------------

final <- create_sample_df(group = "both", tp = "all")

# calculate summary of QA -------------------------------------------------

for (i in 1:nrow(final)) {
  tmp = strsplit(toString(final[i, "scan_dir"]), '/')[[1]][8]
  #if (!stringr::str_detect(tmp, "^ADI\\d{3}_(bl|fu|fu2)$")){
  #  next
  #}
  conf = read.csv(
    paste0("/data/pt_02161/preprocessed/resting/detailedQA/metrics/",
           tmp, "/confounds.csv"))
  tmp = cor.test(conf$DVARS, conf$FD)
  final[i, "corr_DVARS_FD"] = tmp$estimate
  final[i, "mean_DVARS"] = mean(conf$DVARS, na.rm = TRUE)
  final[i, "max_DVARS"] = max(conf$DVARS, na.rm = TRUE)
  final
}
saveRDS(object = final, file = "../report/tab/final.rds")


saveRDS(object = final, file = "../report/tab/final.rds")

# indices of first row of each subject
subj_df <- final[match(unique(final$subj.ID), final$subj.ID),] # df with unique subjects
subj_df$Sex <- plyr::mapvalues(subj_df$Sex, from = c("-1", "1"), to = c("f", "m"))
saveRDS(object = subj_df, file = "../report/tab/subj_df.rds")

# completed time points
freqtable <- as.data.frame(table(final$subj.ID))
tpfreq <- table(freqtable$Freq)

# distribution of condition and sex
cs <- addmargins(table(subj_df$Sex,subj_df$condition))

dfw <- tidyr::pivot_wider(data = final,
                          id_cols = c("subj.ID","condition","Sex","Age_BL"),
                          names_from = "tp",
                          values_from = c("BMI","meanFD"))
saveRDS(object = dfw, file = "../report/tab/dfw.rds")


# tableSample -------------------------------------------------------------

tableSample <- mk_SampleTable(final)
saveRDS(object = tableSample, file = "../report/tab/tableSample.rds")


# tableBaselineCharacteristics --------------------------------------------

dfBL <- create_sample_df(group = "all", tp = "BL")
saveRDS(object = dfBL, file = "../report/tab/dfBL.rds")


# tableDescr --------------------------------------------------------------

# table for descriptive statistics
tableDescr <- psych::describeBy(dfw[,5:ncol(dfw)], group=dfw$condition,mat=T)
tableDescr$predictor <- c(rep("BMI",6),rep("meanFD",6))
tableDescr$tp <- car::recode(tableDescr$vars, "c(1,4)='fu2';c(2,5)='fu';c(3,6)='bl'")
tableDescr$group <- tableDescr$group1
tableDescr <- tableDescr[,c("predictor","tp","group","n","mean","sd")]
tableDescr <- tableDescr[with(tableDescr, order(predictor, tp)),]
tableDescr[,c("mean","sd")] <- round(tableDescr[,c("mean","sd")],2)
tableDescr[c(2:6,8:12),"predictor"] <- ""
tableDescr[c(2,4,6,8,10,12),"tp"] <- ""
rownames(tableDescr) <- NULL
saveRDS(object = tableDescr, file = "../report/tab/tableDescr.rds")


# figBMIdescr -------------------------------------------------------------

fig_BMIdescr <- mk_figBMIdescr(final)
saveRDS(object = fig_BMIdescr, file = "../report/fig/fig_BMIdescr.rds")


# figDesignmatrix ---------------------------------------------------------

DesignMatricesList <- mk_figDesignMatrix()
saveRDS(object = DesignMatricesList, file = "../report/fig/DesignMatricesList.rds")


# mk_tableQcfc ------------------------------------------------------------

res <- read.csv("../report/FDFCcorrelations.csv")
rownames(res) <- c("minimally processed", "AROMA", "AROMA + CC", "AROMA + CC + GSR")
colnames(res) <- c('mean FD-QC', 'median FD-QC', 'sig. vertex', 'sig. vertex BH', 'distance-FD-QC', 'pval')

tableQcfc <- mk_tableQcfc(res)
saveRDS(object = tableQcfc, file = "../report/fig/tableQcfc.rds")

