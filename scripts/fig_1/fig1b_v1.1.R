#The purpose of this script is to create Fig 1b dotplot

library(tidyverse)
library(ggthemes)
library(scales)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(stats)
library(plyr)
library("writexl")
library(PMCMRplus)
library("ggbeeswarm")
library(rstatix)

#import dataset
library(readxl)
#cellseed_ALL <- read_excel("Documents/Seeding Density data/190131-KLG-120-60x-Jurkat-BJAB-THP1-celldensity_20190131_131943[3530]/AssayPlate_PerkinElmer_CellCarrier-384 Ultra[3680]/2019-01-31T131944-0500[3730]/190131-EXP020-53BP1-Spots(260277).result.1.xlsx")


cellseed_ALL <- read_excel("data/raw-data_v1/fig_1/fig_1c/190131-EXP020-53BP1-Spots(260277).result.1.xlsx")
View(cellseed_ALL)
#-------------------------------------------------------------------------
#Process raw dataset to label cell type and seeding density 


#Filter out data that is not needed for analysis 
cellseed_fil <- cellseed_ALL %>% 
  filter(Row %in% c(9:12))

View(cellseed_fil)

cellseed_lab <- cellseed_fil %>% 
  #FIRST ASSIGN DONOR 
  mutate(cell_type = case_when(
    #Donor1
    Column %in% c(3:7) ~"Jurkat"
    #Donor2
    ,Column %in% c(10:14) ~"BJAB"
    #Donor3
    ,Column %in% c(17:21) ~"THP-1"),
    seed_dense = case_when(
      Column %in% c(3,10,17) ~ "20000", 
      Column %in% c(4,11,18) ~"40000", 
      Column %in% c(5,12,19) ~ "60000", 
      Column %in% c(6,13,20) ~ "80000", 
      Column %in% c(7,14,21) ~ "100000")
  )

View(cellseed_lab)

#Set order of seeding density and cell type
cellseed_lab$seed_dense <- factor(cellseed_lab$seed_dense, levels = c("20000","40000","60000","80000","100000"))
cellseed_lab$cell_type <- factor(cellseed_lab$cell_type, levels = c("Jurkat","BJAB", "THP-1"))
#------------------------------------
#SUMMARIZE DATA 
#prepare data for jurkat 
jur_seed <- filter(cellseed_lab, cell_type == "Jurkat")
View(jur_seed)

labels_cell <-cellseed_lab %>% 
  group_by(cell_type,seed_dense) %>%
  dplyr::summarise(counts = n())

#Find range of final nuclei for each cell type 
aggregate(nuc_total~WellName, cellseed_lab, mean)
#%>% 
#summarise(N=paste0("n =", n())) 


#summarise(counts = n()) %>%
View(cellseed_lab)

#-----------------------------------
#Create dotplot
#prepare boxplot w/ dots 
all_dotplot <- ggplot(cellseed_lab, aes(x = seed_dense, y = nuc_total, fill = seed_dense)) +
  geom_dotplot(binaxis='y', binwidth = 15, stackdir='center', dotsize=2.5, show.legend = FALSE) +
  geom_point(stat="summary", fun.y="mean", show.legend = FALSE) + 
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0) +
  #geom_text(data = labels_cell,aes(x=seed_dense, y=2500,label=N), size=6) + 
  ylab("Number of recovered cells") + 
  xlab("Seeding Density"~"(x10"^"4"~"cells/well)") + 
  labs(fill="Seeding Density ") + 
  scale_fill_brewer(palette = "BuGn") + 
  theme_bw(base_size = 20) 

all_dotplot + facet_wrap(~cell_type) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=15), panel.border = element_blank(), panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.background = element_blank(), axis.title.x = element_blank())  + scale_x_discrete(labels = c("20000"= "2", "40000"= "4", "60000"="6", "80000"= "8", "100000" = "10"))



all_dotplot + facet_wrap(~cell_type) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=15), panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(), axis.line = element_blank(), strip.background = element_blank())  + scale_x_discrete(labels = c("20000"= "2", "40000"= "4", "60000"="6", "80000"= "8", "100000" = "10"))


theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank())
# stat_compare_means(comparisons = my.comparisons, label = "p.adj.signif", method = "t.test", paired = FALSE) +    
# annotate("text", x=2.8, y=2, label = "ANOVA, p = 0.99") #t.test comparison 


