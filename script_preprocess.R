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
library('readr')

load('rawpcr.RData')

DeltaCT <- function(ct) {
  ## INPUT: `ct` is numeric tbl_df. The first row is house keeping gene. The first column is the marker annotation.
  ## USAGE: Calculate the delta_ct.
  ## OUTPUT: A tbl_df

  require('magrittr')
  require('dplyr')

  res <- ct %>%
    mutate_at(-1, ~funs(. - first(.))) %>%
    slice(-1)

  return(res)
}

DeltaDeltaCT <- function(deltact) {
  ## INPUT: `deltact` is numeric tbl_df. The last rows are controls. The 1st column the Time, the 2nd column is the marker annotation.
  ## USAGE: Calculate the delta_ct.
  ## OUTPUT: A tbl_df

  require('magrittr')
  require('dplyr')

  res <- deltact %>%
    mutate_at(-1:-2, ~funs(-(. - .[9:12]))) %>%
    slice(-9:-12)

  return(res)
}

GroupMean <- function(ddct) {
  ## INPUT: `ddct` is the deltadeltaCt. The 1st column is label. P1/P3/P5/P7/P9/P11 are control, and P2/P4/P6/P8/P10 are treatment.

  require('magrittr')
  require('dplyr')
  require('tibble')

  control <- ddct %>%
    select(c(2, 4, 6, 8, 10, 12)) %>%
    rowMeans

  treatment <- ddct %>%
    select(c(3, 5, 7, 9, 11)) %>%
    rowMeans

  meanM <- tibble(control = control,
                  treatment = treatment)

  res <- ddct %>%
    select_at(1) %>%
    bind_cols(meanM)

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
    kmeans(., centers = 4, nstart = 25) %>%
    fviz_cluster(., data = plotMat) +
    ggtitle(p[i])

}

ml <- marrangeGrob(resp, nrow = 2, ncol = 1)
ggsave('deltaCt.pdf', ml, width = 21, height = 12)

write.csv(deltaCt, file = 'deltaCt.csv', row.names = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~delta delta Ct~~~~~~~~~~~~~~~~~~~~
deltadeltaCt <- deltaCt %>%
  group_by(., Location) %>%
  do(., DeltaDeltaCT(select(., Time:P11))) %>%
  ungroup

## kmeans
p <- paste0('P', 1:11)
resp <- vector('list', 11)
for (i in seq_along(p)) {

  plotMat <- deltadeltaCt %>%
    select(., Location:Bacteria, p[i]) %>%
    spread(., Location, p[i]) %>%
    mutate(., Label = paste(Bacteria, Time, sep = '_')) %>%
    select(., -Time, -Bacteria) %>%
    as.data.frame %>%
    column_to_rownames(., var = 'Label')

  resp[[i]] <- plotMat %>%
    kmeans(., centers = 4, nstart = 25) %>%
    fviz_cluster(., data = plotMat) +
    ggtitle(p[i])

}

ml <- marrangeGrob(resp, nrow = 2, ncol = 1)
ggsave('negdeltadeltaCt.pdf', ml, width = 21, height = 12)

write.csv(deltadeltaCt, file = 'negdeltadeltaCt.csv', row.names = FALSE)

## dot plot
tmp1 <- deltadeltaCt %>%
  select(., Location:Bacteria, P1) %>%
  spread(., Time, P1) %>%
  ggplot(., aes(x = `3weeks`, y = `3months`)) +
  geom_point(aes(colour = Bacteria), size = 3)

## merge all dot plot
deltadeltaCt %>%
  gather(., P1:P11, key = 'Person', value = 'ddCt') %>%
  spread(., Time, ddCt) %>%
  ggplot(., aes(x = `3weeks`, y = `3months`)) +
  geom_point(aes(colour = Bacteria, shape = Location), size = 3)
ggsave('weekmonthmerge.pdf')


## merge split dot plot
deltadeltaCt %>%
  gather(., P1:P11, key = 'Person', value = 'ddCt') %>%
  spread(., Time, ddCt) %>%
  ggplot(., aes(x = `3weeks`, y = `3months`)) +
  geom_point(aes(colour = Bacteria, shape = Location), size = 3) +
  facet_wrap(. ~ Person, ncol = 2)
ggsave('weekmonthsplit.pdf')

##~~~~~~~~~~~~~~~~~~~~plot and test deltaCt~~~~~~~~~~~~~~~~~~~~~~
ctPlotData <- deltadeltaCt %>%
  gather(P1 : P11, key = 'Sample', value = 'nddct') %>%
  mutate(Group = rep(c(rep(c('control', 'treatment'), 5), 'control'),  each = 32))

ctPlotData %>%
  ggplot(aes(x = Time, y = nddct, fill = Group, colour = Group)) +
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_violin(alpha = 0.3, position = position_dodge(width=0.9)) +
  geom_boxplot(alpha = 0.6, width = 0.1, outlier.colour = NA, position = position_dodge(width=0.9)) +
  facet_wrap(. ~ Bacteria, ncol = 2) +
  xlab('') +
  ylab(expression(paste('-', Delta, Delta, 'Ct')))
ggsave('control_treat_split.pdf')

ctPlotData %>%
  ggplot(aes(x = Time, y = nddct, fill = Group, colour = Group)) +
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_violin(alpha = 0.3, position = position_dodge(width=0.9)) +
  geom_boxplot(alpha = 0.6, width = 0.1, outlier.colour = NA, position = position_dodge(width=0.9)) +
  xlab('') +
  ylab(expression(paste('-', Delta, Delta, 'Ct')))
ggsave('control_treat.pdf')

ctPlotData %>%
  group_by(Time, Group) %>%
  transmute(t_test = nddct %>% t.test %>% .$p.value,
            w_test = nddct %>% wilcox.test %>% .$p.value) %>%
  distinct

ctPlotData %>%
  group_by(Time, Group, Bacteria) %>%
  transmute(t_test = nddct %>% t.test %>% .$p.value,
            w_test = nddct %>% wilcox.test %>% .$p.value) %>%
  distinct

ctPlotData %>%
  group_by(Time) %>%
  transmute(t_test = split(nddct, Group) %>% {t.test(.$control, .$treatment)} %>%.$p.value,
            w_test = split(nddct, Group) %>% {wilcox.test(.$control, .$treatment)} %>%.$p.value) %>%
  distinct

ctPlotData %>%
  group_by(Time, Bacteria) %>%
  transmute(t_test = split(nddct, Group) %>% {t.test(.$control, .$treatment)} %>%.$p.value,
            w_test = split(nddct, Group) %>% {wilcox.test(.$control, .$treatment)} %>%.$p.value) %>%
  distinct
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################################################
