library(kableExtra)

# processes data to report later in manuscript .Rmd file
ROOT_DIR = "/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /code_and_manuscript/code/"

# Get functions -----------------------------------------------------------
setwd(ROOT_DIR)
source("figures.R")
source("tables.R")
source("create_sample_df.R")
source("../../../../../Preprocessing/qa/rs_qa/group_level_QA/QC_FC_correlations/calc_qc_fc_correlations.R")

options(knitr.table.format = "latex")
options(knitr.kable.NA = ' ') # hide NA in tables

# Load data ---------------------------------------------------------------

final <- create_sample_df(group = "both", tp = "all")


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
#write.csv(file = "../report/all_QA_measures.csv", final)

final <- within(final,  subj.ID_tp <- paste(subj.ID,tp, sep="_"))

# save final data frame
write.csv(final,"../report/final.csv", row.names = FALSE)

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
saveRDS(object = dfBL, file = "../report/dfBL.rds")


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

# FDFC correlation computation --------------------------------------------

# ############################ IMPORTANT ##################################
# ONLY COMPUTE IF NECESSARY (long computation time)
# #########################################################################

#FDFC=calc_qc_fc(final)
#write.csv(FDFC,"../report/FDFCcorrelations.csv", row.names = FALSE)


# prepare tsnr ROI plot ---------------------------------------------------

fig_tSNR <- mk_figtSNR(final)
ggsave("../report/tsnr.pdf", height = 6.604, units="cm", fig_tSNR)
#saveRDS(object = fig_tSNR, file = "../report/fig/fig_tSNR.rds")


# load and prepare aggFC measures -----------------------------------------

agg_FC_DMN=read.table("/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/DMN/agg_FC_DMN.txt")
agg_FC_DMN_ID=read.table("/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/DMN/agg_FC_DMN_ID.txt")
agg_FC_Rew_ID=read.table("/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/Rew/agg_FC_Rew_ID.txt")
agg_FC_Rew=read.table("/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis/aggFC/Rew/agg_FC_Rew.txt")
for (i in c(1:nrow(agg_FC_DMN_ID))){
  agg_FC_Rew[i,"subj.ID_tp"]=strsplit(base::strsplit(as.character(agg_FC_Rew_ID[i,"V1"]),'/')[[1]][10],'_nacc')[[1]][1]
  agg_FC_DMN[i,"subj.ID_tp"]=strsplit(base::strsplit(as.character(agg_FC_DMN_ID[i,"V1"]),'/')[[1]][10],'_pcc')[[1]][1]
}

colnames(agg_FC_Rew)=c("mean_Rew_conn","sd_Rew_conn","subj.ID_tp")
colnames(agg_FC_DMN)=c("mean_DMN_conn","sd_DMN_conn","subj.ID_tp")

final_FC=merge(final, agg_FC_Rew, by="subj.ID_tp", all.x=TRUE)
final_FC=merge(final_FC, agg_FC_DMN, by="subj.ID_tp", all.x=TRUE)
final_FC$subj.ID=as.factor(final_FC$subj.ID)
final_FC$group_factor=relevel(final_FC$group_factor, ref = "KG")

# alculate mean BMI over timepoints
final_FC.meanBMI=tapply(X=final_FC$BMI,
                        INDEX=final_FC$subj.ID, FUN=mean, na.rm=TRUE)
final_FC$mean.BMI=
  final_FC.meanBMI[as.numeric(final_FC$subj.ID)]
final_FC$within.BMI=
  final_FC$BMI-final_FC$mean.BMI

# calculate mean FD over timepoints
final_FC.meanFD=tapply(X=final_FC$logmFD,
                       INDEX=final_FC$subj.ID, FUN=mean, na.rm=TRUE)
final_FC$mean.logmFD=
  final_FC.meanFD[as.numeric(final_FC$subj.ID)]
final_FC$within.logmFD=
  final_FC$logmFD-final_FC$mean.logmFD

write.csv(final_FC,"../report/final_FC.csv", row.names = FALSE)


# figDvarsmFD -------------------------------------------------------------

figDvarsmFDList <- figDvarsmFD(final)
saveRDS(object = figDvarsmFDList, file = "../report/fig/figDvarsmFDList.rds")


# fig_DMNconn -------------------------------------------------------------

# DMN conn over three time points
fig_DMNconn <- mk_figDMNdescr(final_FC)
saveRDS(object = fig_DMNconn, file = "../report/fig/figDMNconn.rds")

# fig_RNconn --------------------------------------------------------------

# RN conn over three time points
fig_Rewconn <- mk_figRewdescr(final_FC)
saveRDS(object = fig_Rewconn, file = "../report/fig/figRewconn.rds")


# FCTables ----------------------------------------------------------------

FCTableList <- mk_FCTables()
saveRDS(object = FCTableList, file = "../report/tab/FCTableList.rds")

