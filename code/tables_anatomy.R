
library(dplyr)
DIR_ANALYSIS <- "/data/pt_02161/Results/Project2_resting_state/connectivity/Analysis"
df <- file.path(DIR_ANALYSIS,"noExclFD/result_report.txt")
AnatomyResults <- read.table(df, sep="\t",col.names = c("model","dir","k"),stringsAsFactors = FALSE)


# Construct table title ---------------------------------------------------

AnatomyResults <- tidyr::separate(data = AnatomyResults, col = "model", into = c("model","roi","variable","association"), sep = "_") 

AnatomyResults$start <- "Anatomical labels for clusters functionally connected to the "
AnatomyResults <- 
  tidyr::extract(data = AnatomyResults, col = "roi", into = "denoising", regex = "([a-z]+)", remove = FALSE) %>%
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
  tidyr::extract(data = AnatomyResults, col = "model", into = "adjust", regex = "([a-z]+)", remove = FALSE) %>%
  tidyr::extract(col = "model", into = "model", regex = "([A-Z]+)")

AnatomyResults$model <- AnatomyResults$model %>%
  stringr::str_replace_all('BMIFD', 'BMI-FD')

AnatomyResults$adjust <- AnatomyResults$adjust %>%
  stringr::str_replace_all('agesexfd', 'age, sex and mFD') %>%
  stringr::str_replace_all('agesex', 'age and sex')

AnatomyResults$roi <- AnatomyResults$roi %>%
  stringr::str_replace_all('NACC', 'NAcc')

TabLabelTitles <- with(
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

TabLabelTitles[1]


# table design ------------------------------------------------------------

txtList <- list()

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
  name_model <- AnatomyResults$model[r]
  txtList[[name_model]] <- txt
}

txt_colname <- kableExtra::linebreak(c("Cluster","Number of voxels\n in cluster","\\% of cluster volume\nassigned", "Hemisphere","Area","\\% of area overlap\nwith cluster"))

TableList <- list()
for (i in 1:length(txtList)){
  TableList[[i]] <-
    txtList[[i]] %>%
    knitr::kable(
      escape = FALSE,     # use font spec of latex with kableExtra
      col.names = txt_colname,
      row.names = FALSE,
      format = "latex",
      booktabs = T,
      linesep = "",      # disable "\\addlinespace" at every 5th line
      #label = name_model[i],
      caption = TabLabelTitles[i],#sprintf("Anatomical labelling for %s",name_model[i]),
    ) %>%
    kableExtra::kable_styling(latex_options = c('scale_down'))
}


TableList[[1]]
TableList[[2]]
TableList[[3]]
TableList[[4]]
TableList[[5]]
TableList[[6]]
TableList[[7]]
TableList[[8]] 
TableList[[9]]
