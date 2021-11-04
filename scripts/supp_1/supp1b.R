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
X190810_203 <- read_excel("Desktop/exp203_CD4_CD8_B_Monocytes_KG/data_203/190810-203-ALL-revised.xlsx")

#Filter data to DMSO and ETP-30µm treatment groups 
#CD4
m_CD4 <- filter(X190810_203, cell_type == "CD4") 
m_CD4$treat <- factor(m_CD4$treat, levels = c("Untreated","ETP"))
m_CD4_perc <- m_CD4 %>%
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#CD8
m_CD8 <- filter(X190810_203, cell_type == "CD8") 
#Order two treatment groups 
m_CD8$treat <- factor(m_CD8$treat, levels = c("Untreated","ETP"))
#percentage filter 
m_CD8_perc <- m_CD8 %>%
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#B-cell
m_B <- filter(X190810_203, cell_type == "B-cell")
m_B$treat <- factor(m_B$treat, levels = c("Untreated","ETP"))
High_spots_yH2AX_ <- filter(m_B,treat_yH2AX == "High")
#percentage 
m_B_perc <- m_B %>% 
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))


#CD4---------------------------------
#Analyze 53BP1 data as a percentage 
j_fil_Jurperc <- j_fil_Jur %>%
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#Create barplot 
bar_CD4 <- ggplot(m_CD4_perc, aes(x = treat, y = perc*100)) +
  geom_bar(
    aes(fill = factor(treat_lev, levels=c("Low","High"))),
    stat = "identity", position = position_stack(), show.legend = TRUE
  ) + 
  ggtitle("CD4+ T cells") + 
  labs(x = "Treatment", y = "Percentage of Total Cells (%)", fill = "53BP1 Spot Status") 
#geom_text(
#aes(label = perc*100, group = treat_lev)

#print barplot
CD4_bar <- bar_CD4 + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text = element_text(size = 15), axis.title.x = element_blank())


#CD8-------------------------------

#Create barplot 
bar_CD8 <- ggplot(m_CD8_perc, aes(x = treat, y = perc*100)) +
  geom_bar(
    aes(fill = factor(treat_lev, levels=c("Low","High"))),
    stat = "identity", position = position_stack(), show.legend = TRUE
  ) + 
  ggtitle("CD8+ T cells") + 
  labs(x = "Treatment", y = "Percentage of Total Cells (%)", fill = "53BP1 Spot Status") 

#print barplot 
CD8_bar <- bar_CD8 + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text = element_text(size = 15), axis.title.x = element_blank())

#B-cell-------------------------------

#Create barplot 
bar_B <- ggplot(m_B_perc, aes(x = treat, y = perc*100)) +
  geom_bar(
    aes(fill = factor(treat_lev, levels=c("Low","High"))),
    stat = "identity", position = position_stack(), show.legend = TRUE
  ) + 
  ggtitle("B cells") + 
  labs(x = "Treatment", y = "Percentage of Total Cells (%)", fill = "53BP1 Spot Status") 

#print barplot 
B_bar <- bar_B + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text = element_text(size = 15), axis.title.x = element_blank())


#Create figure-------------

supp1b <- ggarrange(CD4_bar, CD8_bar, B_bar, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom")

supp1b
