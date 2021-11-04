#The purpose of this script is to visualize supplement 1a 

library(tidyverse)
library(ggthemes)
library(scales)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(stats)
library(forcats)
library(ggsignif)


#find and read data from the directory
library(readxl)
setwd("/Users/gallantkl")

#Load data 
X_122_jurkat <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190218-121-KG-jurkat-yH2AX-53bp1-dose_20190218_141455[3607]/190221-EXP122-DoseResponse-Jurkat.xlsx")
X_119_BJAB <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190124-119-KLG-Jurkat-BJAB-THP1-Ad2x_20190124_113931[3510]/190218-119-BJAB.xlsx")
X_122_THP1 <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190218-122-KG-thp1-yH2AX-53bp1-dose_20190218_150512[3608]/190218-EXP122-DoseResponse-THP1.xlsx")

#Filter data to DMSO and ETP-30µm treatment groups 
j_fil_Jur <- X_122_jurkat %>%
  filter(treat %in% c("DMSO", "30")) 

j_fil_THP <- X_122_THP1 %>%
  filter(treat %in% c("DMSO", "30"))

j_fil_BJAB <- X_119_BJAB %>%
  filter(treat %in% c("DMSO", "30"))

j_fil_Jur$treat <- factor(j_fil_Jur$treat, levels = c("DMSO","30"))
j_fil_BJAB$treat <- factor(j_fil_BJAB$treat, levels = c("DMSO","30"))
j_fil_THP$treat <- factor(j_fil_THP$treat, levels = c("DMSO","30"))

#Jurkat---------------------------------
#Analyze 53BP1 data as a percentage 
j_fil_Jurperc <- j_fil_Jur %>%
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#Create barplot 
j_barplot <- ggplot(j_fil_Jurperc, aes(x = treat, y = perc*100)) +
  geom_bar(
    aes(fill = treat_lev),
    stat = "identity", position = position_stack(), show.legend = TRUE
  ) + 
  ggtitle("Jurkat") + 
  labs(x = "Treatment", y = "Percentage of Total Cells (%)", fill = "53BP1 Spot Status", main = "Jurkat")

#print barplot
jurkat_bar <- j_barplot +  scale_fill_manual(values=c("#e5f5f8", "#99d8c8"), labels = c('Low', 'High'))  + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text = element_text(size = 15), axis.title.x = element_blank()) 


#BJAB-------------------------------
#Analyze 53BP1 data as a percentage 
j_fil_BJABperc <- j_fil_BJAB %>%
  group_by(treat, treat_lev) %>%
  dplyr::summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#Create barplot 
b_barplot <- ggplot(j_fil_BJABperc, aes(x = treat, y = perc*100)) +
  geom_bar(
    aes(fill = treat_lev),
    stat = "identity", position = position_stack(), show.legend = TRUE
  ) + 
  ggtitle("BJAB") + 
  labs(x = "Treatment", y = "Percentage of Total Cells (%)", fill = "53BP1 Spot Status")  

#print barplot 
bjab_bar <- b_barplot + scale_fill_manual(values=c("#e5f5f8", "#99d8c8"), labels = c('Low', 'High'))  + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text = element_text(size = 15), axis.title.x = element_blank()) 


#THP1--------------------------------
#Analyze 53BP1 data as a percentage 
j_fil_THPperc <- j_fil_THP %>%
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#Create barplot 
t_bar <- ggplot(j_fil_THPperc, aes(x = treat, y = perc*100)) +
  geom_bar(aes(fill = treat_lev),
           stat = "identity", position = position_stack(reverse = TRUE), show.legend = TRUE) + 
  ggtitle("THP-1") +
  labs(x = "Treatment", y = "Percentage of Total Cells (%)", fill = "53BP1 Spot Status") 

#print barplot 
thp1_bar <- t_bar + scale_fill_manual(values=c("#99d8c8","#e5f5f8"), labels = c('Low', 'High'))  + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text = element_text(size = 15), axis.title.x = element_blank()) 


#Create figure-------------

supp1a <- ggarrange(jurkat_bar, bjab_bar, thp1_bar, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom")

supp1a
