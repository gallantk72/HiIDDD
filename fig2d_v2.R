#The purpose of this script is to plot figure 2d 
#load libraries 
#the purpose of this code is for data visualization of cells 
library(tidyverse)
library(ggthemes)
library(scales)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(stats)
library(forcats)
library(ggsignif)
library(readxl)

#set working directory 
setwd("/Users/gallantkl")

#load datasets 
X_122_jurkat <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190218-121-KG-jurkat-yH2AX-53bp1-dose_20190218_141455[3607]/190221-EXP122-DoseResponse-Jurkat.xlsx")
X_119_BJAB <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190124-119-KLG-Jurkat-BJAB-THP1-Ad2x_20190124_113931[3510]/190218-119-BJAB.xlsx")
X_122_THP1 <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190218-122-KG-thp1-yH2AX-53bp1-dose_20190218_150512[3608]/190218-EXP122-DoseResponse-THP1.xlsx")

#prepare datasets 
#filter dataset to specific density 
j_fil_Jur <- X_122_jurkat %>%
  filter(treat %in% c("DMSO", "30")) 

j_fil_BJAB <- X_119_BJAB %>%
  filter(treat %in% c("DMSO", "30"))

View(j_fil_BJAB)


j_fil_THP <- X_122_THP1 %>%
  filter(treat %in% c("DMSO", "30"))


#plot confilled 


#plot jurkat------
Jurkat_confilled <- ggplot(j_fil_Jur, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat), show.legend = FALSE) +
  xlab("53BP1 Sum Spot Intensity (A.U.)") + 
  ylab(expression(gamma ~ "H2AX Mean Nuclear Intensity (A.U.)")) + 
  labs(fill="Treatment") + 
  ggtitle("Jurkat") + 
  scale_x_continuous(trans = 'log2', labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific)


fig2d_jurkat <- Jurkat_confilled + theme_bw(base_size = 20) +
  theme(axis.text = element_text(size=15), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP"))

#plot bjab ------
BJAB_confilled <- ggplot(j_fil_BJAB, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat), show.legend = FALSE) +
  xlab("53BP1 Sum Spot Intensity (A.U.)") + 
  ylab(expression(gamma ~ "H2AX Mean Nuclear Intensity (A.U.)")) + 
  labs(fill="Treatment") + 
  ggtitle("BJAB") + 
  scale_x_continuous(trans = 'log2', labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific)


fig2d_bjab <- BJAB_confilled + theme_bw(base_size = 20) +
  theme(axis.text = element_text(size=15), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP"))

#plot thp1 -----
THP1_confilled <- ggplot(j_fil_THP, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat), show.legend = FALSE) +
  xlab("53BP1 Sum Spot Intensity (A.U.)") + 
  ylab(expression(gamma ~ "H2AX Mean Nuclear Intensity (A.U.)")) + 
  labs(fill="Treatment") + 
  ggtitle("THP-1") + 
  scale_x_continuous(trans = 'log2', labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific)


fig2d_thp1 <- THP1_confilled + theme_bw(base_size = 20) +
  theme(axis.text = element_text(size=15), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP"))


#plot fig2d ---------
fig2d <- ggarrange(fig2d_jurkat, fig2d_bjab, fig2d_thp1,
                   labels = c("a", "b", "c"),  
                   ncol = 3, nrow = 1)
fig2d
