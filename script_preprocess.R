######################deal with NA#################
setwd('/extDisk2/toothepaste/')

library('readr')
library('magrittr')
library('dplyr')

rawpcr <- read_csv('normalize.csv')

## use average in group
## rawpcr$P1[8] <- mean(rawpcr$P1[7:9], na.rm=TRUE)
## rawpcr$P2[124] <- mean(rawpcr$P2[124:126], na.rm=TRUE)
## rawpcr$P3[95] <- mean(rawpcr$P3[94:96], na.rm=TRUE)
## rawpcr$P4[111] <- mean(rawpcr$P4[109:111], na.rm=TRUE)
## rawpcr$P4[123] <- mean(rawpcr$P4[121:123], na.rm=TRUE)
## rawpcr$P4[143] <- mean(rawpcr$P4[142:144], na.rm=TRUE)
## rawpcr$P5[66] <- mean(rawpcr$P5[64:66], na.rm = TRUE)
## rawpcr$P5[68] <- mean(rawpcr$P5[67:69], na.rm = TRUE)
## rawpcr$P6[33] <- mean(rawpcr$P6[31:33], na.rm = TRUE)
## rawpcr$P9[7] <- mean(rawpcr$P9[7:9], na.rm = TRUE)
## rawpcr$P9[17:18] <- mean(rawpcr$P9[16:18], na.rm = TRUE)
## rawpcr$P9[60] <- mean(rawpcr$P9[58:60], na.rm = TRUE)
## rawpcr$P9[121] <- mean(rawpcr$P9[121:123], na.rm = TRUE)
## rawpcr$P10[8] <- mean(rawpcr$P10[7:9], na.rm = TRUE)
## rawpcr$P10[41] <- mean(rawpcr$P10[40:42], na.rm = TRUE)
## rawpcr$P10[68] <- mean(rawpcr$P10[67:69], na.rm = TRUE)
## rawpcr$P11[18] <- mean(rawpcr$P11[16:18], na.rm = TRUE)
## rawpcr$P11[22] <- mean(rawpcr$P11[22:24], na.rm = TRUE)
## rawpcr$P11[38] <- mean(rawpcr$P11[37:39], na.rm = TRUE)
## rawpcr$P11[67] <- mean(rawpcr$P11[67:69], na.rm = TRUE)
## rawpcr$P11[98] <- mean(rawpcr$P11[97:99], na.rm = TRUE)
## rawpcr$P11[127] <- mean(rawpcr$P11[127:129], na.rm = TRUE)

## blank group use delta-CT
ctAA3W <- rawpcr %>%
  filter(., Location == 'A', Time == '3weeks', Bacteria %in% c('AA', '16S')) %>%
  select(., -P2, -Location, -Time) %>%
  group_by(., Bacteria) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  select(., -Bacteria)

rawpcr$P2[64:66] <- ((slice(ctAA3W, 2) - slice(ctAA3W, 1)) %>% unlist %>% mean) + mean(rawpcr$P2[61:63])

save(rawpcr, file = 'rawpcr.RData')
################################################


##########################manipulate CT#####################
setwd('/extDisk2/toothepaste/')

library('dplyr')
library('tidyr')
library('magrittr')
library('ggplot2')
library('tibble')
library('factoextra')
library('gridExtra')

load('rawpcr.RData')

DeltaCT <- function(ct) {
  ## INPUT: CT is numeric tbl_df. The first row is house keeping gene. The first column is the marker annotation.
  ## USAGE: Calculate the delta_ct.
  ## OUTPUT: A tbl_df

  require('magrittr')
  require('dplyr')

  res <- ct %>%
    mutate_at(., -1, funs(. - first(.))) %>%
    slice(-1)

  return(res)
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~delta Ct~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
deltaCt <- rawpcr %>%
  group_by(., Location, Time, Bacteria) %>%
  summarise_all(., mean, na.rm = TRUE) %>%
  do(., DeltaCT(select(., Bacteria:P11))) %>%
  ungroup

##  A B C D
## AA_3months
## AA_3weeks
p <- paste0('P', 1:11)
resp <- vector('list', 11)
for (i in seq_along(p)) {

  plotMat <- deltaCt %>%
    select(., Location:Bacteria, p[i]) %>%
    spread(., Location, p[i]) %>%
    mutate(., Label = paste(Bacteria, Time, sep = '_')) %>%
    select(., -Time, -Bacteria) %>%
    as.data.frame %>%
    column_to_rownames(., var = 'Label')

  resp[[i]] <- plotMat %>%
    kmeans(., centers = 3, nstart = 25) %>%
    fviz_cluster(., data = plotMat) +
    ggtitle(p[i])

}

ml <- marrangeGrob(resp, nrow = 2, ncol = 1)
ggsave('deltaCt.pdf', ml, width = 21, height = 12)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~delta deltaCt~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################################################
