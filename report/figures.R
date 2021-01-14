
pink_blue = c("#00305E", "#C00045")
blue_organge = c("#046C9A", "#D69C4E")

# BMI --------------------------------------------------------------------------
mk_figBMIdescr <- function(final) {

  condition.labs <- c("BARS","NBARS")
  names(condition.labs) <- c("IG","KG")
  # BMI over time each group
  figBMIdescr <-
    ggplot(final, aes(
      x = tp,
      y = BMI,
      group = subj.ID,
      color = condition
    )) +
    geom_line() + geom_point() +
    stat_summary(
      aes(group = condition),
      fun = mean,
      geom  = "line",
      size = 2,
      color = c(rep("#046C9A", 3), rep("#D69C4E", 3))
    ) +
    facet_grid(. ~ condition, labeller = labeller(condition=condition.labs)) +
    xlab("time points") +
    scale_x_discrete(labels = c("bl" = "0", "fu" = "6", "fu2" = "12")) +
    scale_color_manual(
      values = c("#046C9A70", "#D69C4E70"),
      labels = c("intervention", "control")
    ) +
    xlab("Time after intervention [months]") + ylab ("BMI [kg/m²]") +
    theme_bw() +  theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "none"
    )
  # ggsave("/data/pt_02161/Publications/Abstracts/Heinrichs_OHBM2021/bmi.pdf", units = "cm", width=15, height=8)
  return(figBMIdescr)
}

# Reward -----------------------------------------------------------------------
mk_figRewdescr <- function(final_FC) {

  condition.labs <- c("BARS","NBARS")
  names(condition.labs) <- c("IG","KG")
  # BMI over time each group
  figRewdescr <-
    ggplot(final_FC, aes(
      x = tp,
      y = mean_Rew_conn,
      group = subj.ID,
      color = condition
    )) +
    geom_line() + geom_point() +
    stat_summary(
      aes(group = condition),
      fun = mean,
      geom  = "line",
      size = 2,
      color = c(rep("#046C9A", 3), rep("#D69C4E", 3))
    ) +
    facet_grid(. ~ condition, labeller = labeller(condition=condition.labs)) + xlab("time points")
  figRewdescr=figRewdescr +
    scale_x_discrete(labels = c("bl" = "0", "fu" = "6", "fu2" = "12")) +
    scale_color_manual(
      values = c("#046C9A70", "#D69C4E70"),
      labels = c("intervention", "control")
    ) +
    xlab("Time after intervention [months]") + ylab ("mean reward connectivity") +
    theme_bw() +  theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "none"
    )
  return(figRewdescr)
}

# DMN --------------------------------------------------------------------------
mk_figDMNdescr <- function(final_FC) {

  condition.labs <- c("BARS","NBARS")
  names(condition.labs) <- c("IG","KG")
  # BMI over time each group
  figDMNdescr <-
    ggplot(final_FC, aes(
      x = tp,
      y = mean_DMN_conn,
      group = subj.ID,
      color = condition
    )) +
    geom_line() + geom_point() +
    stat_summary(
      aes(group = condition),
      fun = mean,
      geom  = "line",
      size = 2,
      color = c(rep("#046C9A", 3), rep("#D69C4E", 3))
    ) +
    facet_grid(. ~ condition, labeller = labeller(condition=condition.labs)) + xlab("time points")
  figDMNdescr=figDMNdescr +
    scale_x_discrete(labels = c("bl" = "0", "fu" = "6", "fu2" = "12")) +
    scale_color_manual(
      values = c("#046C9A70", "#D69C4E70"),
      labels = c("intervention", "control")
    ) +
    xlab("Time after intervention [months]") + ylab ("mean DMN connectivity") +
    theme_bw() +  theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "none"
    )
  return(figDMNdescr)
}

# Design matrix ----------------------------------------------------------------
mk_figDesignMatrix <- function() {

  # design matrix time as factor
  txt_path <- "/data/pt_02161/Analysis/Project2_resting_state/seed-based/Second_level /SwE_files/noExclFD/both/total/"

  ID <- read.delim(paste0(txt_path,"subjID.txt"), head=FALSE)
  condition <- read.delim(paste0(txt_path,"group.txt"), head=FALSE)
  KG_bl <- read.delim(paste0(txt_path,"KG_bl.txt"), head=FALSE)
  KG_fu <- read.delim(paste0(txt_path,"KG_fu.txt"), head=FALSE)
  KG_fu2 <- read.delim(paste0(txt_path,"KG_fu2.txt"), head=FALSE)
  IG_bl <- read.delim(paste0(txt_path,"IG_bl.txt"), head=FALSE)
  IG_fu <- read.delim(paste0(txt_path,"IG_fu.txt"), head=FALSE)
  IG_fu2 <- read.delim(paste0(txt_path,"IG_fu2.txt"), head=FALSE)
  measurement <- c(1:nrow(ID))
  dmf <- data.frame(measurement,ID,condition,KG_bl,KG_fu,KG_fu2,IG_bl,IG_fu,IG_fu2)
  colnames(dmf) <- c("measurement","ID","condition","N0","N6","N12","B0","B6","B12")

  # design matrix time as continuous variable
  ID <- read.delim(paste0(txt_path,"subjID.txt"), head=FALSE)
  condition <- read.delim(paste0(txt_path,"group.txt"), head=FALSE)
  group_IG <- read.delim(paste0(txt_path,"group_IG.txt"), head=FALSE)
  group_KG <- read.delim(paste0(txt_path,"group_KG.txt"), head=FALSE)
  tp_IG <- read.delim(paste0(txt_path,"tp_IG.txt"), head=FALSE)
  tp_KG <- read.delim(paste0(txt_path,"tp_KG.txt"), head=FALSE)
  measurement <- c(1:nrow(ID))
  dmc <- data.frame(measurement,ID,condition,group_IG,group_KG,tp_IG,tp_KG)
  colnames(dmc) <- c("measurement","ID","condition","groupB", "groupN","timeB","timeN")

  # plot dmf
  input_dmf <- dmf %>%
    pivot_longer(
      cols= colnames(dmf[,c(4:ncol(dmf))]),
      names_to = "regressor")
  input_dmc <- dmc %>%
    pivot_longer(
      cols= colnames(dmc[,c(4:ncol(dmc))]),
      names_to = "regressor")

  plotA <- ggplot(input_dmf, aes(x = regressor, y = measurement, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low="black", high="white") +
    labs(x="regressors", y="scans") +
    # design
    scale_x_discrete(expand = c(0,0),labels=colnames(dmf)[4:ncol(dmf)]) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() + theme(legend.position = "none",
                       axis.text = element_text(size = 10),
                       axis.title = element_text(size = 12),
                       strip.text = element_text(size = 12))
  # + theme(panel.grid = element_blank(),
  #         panel.border = element_rect(colour = "black", size=1),
  #    legend.text = element_text(size = 10)
  #         plot.title=element_text(size=11))
  plotB <- ggplot(input_dmc, aes(x = regressor, y = measurement, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low="black", high="white") +
    labs(x="regressors", y="scans") +
    # design
    scale_x_discrete(expand = c(0,0),labels=colnames(dmc)[4:ncol(dmc)]) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() + theme(legend.position = "none",
                       axis.text = element_text(size = 10),
                       axis.title = element_text(size = 12),
                       strip.text = element_text(size = 12))
  # + theme(panel.grid = element_blank(),
  #         panel.border = element_rect(colour = "black", size=1),
  #    legend.text = element_text(size = 10)
  #         plot.title=element_text(size=11))
  PlotList <- list(plotA, plotB)
  names(PlotList) <- c("plotA", "plotB")
  return(PlotList)
}

mk_figtSNR <- function (final) {
  final$subj.ID_tp=paste0(final$subj.ID,'_', final$tp)
  
  tsnr=read.table("/data/pt_02161/Analysis/Preprocessing/qa/rs_qa/group_level_QA/check_tsnr/tsnr_in_ROIs.txt")
  tsnr=tsnr[tsnr$V1!="subj",]
  colnames(tsnr)=c("subj.ID_tp","roi","median","mean","sd")
  levels(tsnr$roi)=droplevels(tsnr$roi)
  levels(tsnr$roi)=c("NAcc", "precuneus")
  tsnr$mean=as.double(tsnr$mean)
  
  final_tsnr=merge(final, tsnr, by="subj.ID_tp", all.x=TRUE)
  
  p <- ggplot(aes(as.factor(roi), mean),data=final_tsnr) +
    geom_violin(aes(fill=roi)) + geom_jitter(height = 0, width = 0.1) + xlab("Region of interest") + ylab("average tSNR") +
    scale_fill_manual(values=wes_palette("Darjeeling2",4,type="discrete")[2:3]) +
    theme_bw() + theme(axis.text=element_text(size=10),
                       axis.title=element_text(size=12),
                       strip.text = element_text(size=12),
                       legend.position="")
  return(p)
}

# exploratory and example plots ------------------------------------------------

mk_otherplots <- function(final){
  # BMI over time per group
  tspag <-
    ggplot(final, aes(x = tp, y = BMI)) +
    stat_summary(
      aes(group = condition),
      fun = mean,
      geom  = "line",
      size = 2,
      color = c(rep("#046C9A", 3), rep("#D69C4E", 3))
    ) +
    geom_line(aes(group = subj.ID, color = condition)) +
    geom_point(aes(color = condition))
  tspag +
    scale_x_discrete(labels = c("bl" = "0", "fu" = "6", "fu2" = "12")) +
    scale_color_manual(
      values = c("#046C9A70", "#D69C4E70"),
      labels = c("intervention", "control")
    ) +
    xlab("Time after intervention [months]") + ylab ("BMI [kg/m²]") +
    theme_bw() +  theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "bottom",
      legend.background = element_rect(colour = "grey")
    )

  # BMI over time per subject
  # group x time interaction for BMI
  tspag <- ggplot(final, aes(
    x = tp,
    y = BMI,
    color = condition,
    group = condition
  )) +
    geom_point() +

    stat_summary(fun = mean, geom = "point") +
    stat_summary(fun = mean, geom = "line")
  tspag  +
    scale_color_manual(values = c("#046C9A", "#D69C4E")) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.position = "top"
    )

  ## example plots

  df_wide <- tidyr::pivot_wider(
    data = final,
    id_cols = "subj.ID",
    names_from = "tp",
    values_from = c("tp", "BMI", "logmFD")
  )

  if (tp == "BLFU") {
    s = 2
  } else {
    s = 1
  }

  BMI_vector <- c("BMI_fu2", "BMI_fu", "BMI_bl")
  BMI_vector <- c("BMI_bl", "BMI_fu", "BMI_fu2")
  FD_vector <- c("logmFD_fu2", "logmFD_fu", "logmFD_bl")
  df_wide$avgBMI <- apply(df_wide[, BMI_vector[s:3]],
                          MARGIN = 1,
                          FUN = mean,
                          na.rm = TRUE)
  df_wide$avgFD <- apply(df_wide[, FD_vector[s:3]],
                         MARGIN = 1,
                         FUN = mean,
                         na.rm = TRUE)
  df_long <- tidyr::pivot_longer(
    data = df_wide,
    cols = all_of(order_vector),
    # bl, fu, fu2 (maintain original order)
    # CAREFUL: Order according to dataframe df
    names_to = "tp_tp",
    values_to = "tp",
    values_drop_na = TRUE
  )
  all(df_long$tp == final$tp) # check if order correct

  head(final[, c("subj.ID", "tp")])
  head(df_long[, c("subj.ID", "tp")])

  df_long$cgnBMI <- final$BMI - df_long$avgBMI
  df_long$avgBMIc <-
    df_long$avgBMI - mean(df_long$avgBMI, na.rm = TRUE)
  df$avgBMIc <- df_long$avgBMIc
  df$cgnBMI <- df_long$cgnBMI

  df_long$cgnFD <- final$logmFD - df_long$avgFD
  df_long$avgFDc <-
    df_long$avgFD - mean(df_long$avgFD, na.rm = TRUE)
  final$avgFDc <- df_long$avgFDc
  final$cgnFD <- df_long$cgnFD


  # example plots for interaction effect (arbitrary variables for illustration purposes)
  ggplot(final, aes(x=tp, y=BMI, group=condition, color=condition)) +
    #geom_smooth(aes(group=condition)) +
    stat_summary(fun = mean, geom = "point", size=3, shape=15) +
    stat_summary(fun = mean, geom = "line", size=1) +
    #scale_color_brewer(values = "Set1") +
    scale_x_discrete(labels=c("bl"= "0", "fu"= "6", "fu2" = "12")) +
    scale_color_manual(values=c("#046C9A","#D69C4E"),
                       labels = c("intervention", "control")) +
    xlab("Time after intervention [months]") + ylab ("BMI [kg/m²]") +
    theme_bw() +  theme(
      axis.text = element_text(size = 10), axis.title = element_text(size = 12),
      strip.text = element_text(size = 12), legend.text = element_text(size = 10),
      legend.position = "bottom", legend.background = element_rect(colour = "grey"))

  # example plot for change effect (purely descriptive)
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  cgnplot <- ggplot(final, aes(x=a,y=logmFD, group=subj.ID)) +
    geom_line(aes(x=BMI, color=subj.ID))+
    geom_point(aes(x=BMI, color=subj.ID, shape=tp), size=3, alpha = 0.5)+
    facet_grid(. ~ condition)
  cgnplot +
    labs(shape="Time after intervention [months]",labels=c("0","6","12")) +
    guides(color = FALSE, alpha=FALSE, size=FALSE) +
    scale_color_manual(values = getPalette(length(unique(final$subj.ID)))) +
    ylab("log(mean FD) (mm)") + xlab("BMI [kg/m²]") + xlim(30,55) +
    theme_bw() +  theme(
      axis.text = element_text(size = 10), axis.title = element_text(size = 12),
      strip.text = element_text(size = 12), legend.text = element_text(size = 10),
      legend.position = "bottom", legend.background = element_rect(colour = "grey"))

  # example plot for change effect (purely descriptive)
  tspag=ggplot(data=final, aes(y=cgnFD, x=cgnBMI, color=condition))+
    geom_point(size=2)+
    #geom_line(aes(color=subj.ID))+
    #geom_smooth(method = "lm", colour = "black", linetype = 1, fill = "mistyrose3" ) +
    xlab("within-subject BMI change [kg/m²]")+
    ylab("within-subject log mean FD change [mm]")+
    guides(alpha=FALSE, size=FALSE,
           color = guide_legend(override.aes = list(size=4))) +
    facet_grid(. ~ condition)
  tspag + scale_color_manual(values=c("#046C9A","#D69C4E")) +
    theme_bw() + theme(axis.text=element_text(colour = "black",size=10),
                       axis.title=element_text(size=14),
                       axis.text.x = element_text(size=12),
                       axis.text.y = element_text(size=12),
                       strip.text = element_text(size=14),
                       legend.position="top",
                       legend.title = element_text(size = 14),
                       legend.text = element_text(size = 14))


}

figDvarsmFD <- function(final){

  # figures --------------------------------------------------------------------

  condition.labs <- c("BARS","NBARS")
  names(condition.labs) <- c("IG","KG")

  # DVARS over time each group
  fig_DVARSdescr <-
    ggplot(final, aes(
      x = tp,
      y = mean_DVARS,
      group = subj.ID,
      color = condition
    )) +
    geom_line() + geom_point() +
    stat_summary(
      aes(group = condition),
      fun = mean,
      geom  = "line",
      size = 2,
      color = c(rep("#046C9A", 3), rep("#D69C4E", 3))
    ) +
    facet_grid(. ~ condition, labeller = labeller(condition=condition.labs)) + xlab("time points")
  fig1=fig_DVARSdescr +
    scale_x_discrete(labels = c("bl" = "0", "fu" = "6", "fu2" = "12")) +
    scale_color_manual(
      values = c("#046C9A70", "#D69C4E70"),
      labels = c("intervention", "control")
    ) +
    xlab("Time [months]") + ylab ("Mean DVARS") +
    theme_bw() +  theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "none"
    )
  # mFD over time each group
  fig_mFDdescr <-
    ggplot(final, aes(
      x = tp,
      y = meanFD,
      group = subj.ID,
      color = condition
    )) +
    geom_line() + geom_point() +
    stat_summary(
      aes(group = condition),
      fun = mean,
      geom  = "line",
      size = 2,
      color = c(rep("#046C9A", 3), rep("#D69C4E", 3))
    ) +
    facet_grid(. ~ condition, labeller = labeller(condition=condition.labs)) + xlab("time points")
  fig2=fig_mFDdescr +
    scale_x_discrete(labels = c("bl" = "0", "fu" = "6", "fu2" = "12")) +
    scale_color_manual(
      values = c("#046C9A70", "#D69C4E70"),
      labels = c("intervention", "control")
    ) +
    xlab("Time [months]") + ylab ("mFD in mm") +
    theme_bw() +  theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "none"
    )
  # DVARS-FD over time each group
  fig_mFDDVARSdescr <-
    ggplot(final, aes(
      x = tp,
      y = corr_DVARS_FD,
      group = subj.ID,
      color = condition
    )) +
    geom_line() + geom_point() +
    stat_summary(
      aes(group = condition),
      fun = mean,
      geom  = "line",
      size = 2,
      color = c(rep("#046C9A", 3), rep("#D69C4E", 3))
    ) +
    facet_grid(. ~ condition, labeller = labeller(condition=condition.labs)) + xlab("time points")
  fig3=fig_mFDDVARSdescr +
    scale_x_discrete(labels = c("bl" = "0", "fu" = "6", "fu2" = "12")) +
    scale_color_manual(
      values = c("#046C9A70", "#D69C4E70"),
      labels = c("intervention", "control")
    ) +
    xlab("Time [months]") + ylab ("DVARS-FD correlation") +
    theme_bw() +  theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "none"
    )
  FigureList <- list(fig1, fig2,fig3)
  names(FigureList) <- c("fig1","fig2","fig3")

  return(FigureList)
}
