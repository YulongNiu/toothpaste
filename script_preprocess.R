######################deal with NA#################
setwd('/extDisk2/toothepaste/')

library('readr')
library('magrittr')

rawpcr <- read_csv('normalize.csv')

## use average in group
rawpcr$P1[8] <- mean(rawpcr$P1[7:9], na.rm=TRUE)
rawpcr$P2[124] <- mean(rawpcr$P2[124:126], na.rm=TRUE)
rawpcr$P3[95] <- mean(rawpcr$P3[94:96], na.rm=TRUE)
rawpcr$P4[111] <- mean(rawpcr$P4[109:111], na.rm=TRUE)
rawpcr$P4[123] <- mean(rawpcr$P4[121:123], na.rm=TRUE)
rawpcr$P4[143] <- mean(rawpcr$P4[142:144], na.rm=TRUE)
rawpcr$P5[66] <- mean(rawpcr$P5[64:66], na.rm = TRUE)
rawpcr$P5[68] <- mean(rawpcr$P5[67:69], na.rm = TRUE)
rawpcr$P6[33] <- mean(rawpcr$P6[31:33], na.rm = TRUE)
rawpcr$P9[7] <- mean(rawpcr$P9[7:9], na.rm = TRUE)
rawpcr$P9[17:18] <- mean(rawpcr$P9[16:18], na.rm = TRUE)
rawpcr$P9[60] <- mean(rawpcr$P9[58:60], na.rm = TRUE)
rawpcr$P9[121] <- mean(rawpcr$P9[121:123], na.rm = TRUE)
rawpcr$P10[8] <- mean(rawpcr$P10[7:9], na.rm = TRUE)
rawpcr$P10[41] <- mean(rawpcr$P10[40:42], na.rm = TRUE)
rawpcr$P10[68] <- mean(rawpcr$P10[67:69], na.rm = TRUE)
rawpcr$P11[18] <- mean(rawpcr$P11[16:18], na.rm = TRUE)
rawpcr$P11[22] <- mean(rawpcr$P11[22:24], na.rm = TRUE)
rawpcr$P11[38] <- mean(rawpcr$P11[37:39], na.rm = TRUE)
rawpcr$P11[67] <- mean(rawpcr$P11[67:69], na.rm = TRUE)
rawpcr$P11[98] <- mean(rawpcr$P11[97:99], na.rm = TRUE)
rawpcr$P11[127] <- mean(rawpcr$P11[127:129], na.rm = TRUE)

save(rawpcr, file = 'rawpcr.RData')
################################################
