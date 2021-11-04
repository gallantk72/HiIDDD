#The purpose of this code is to generate Fig. 1C from manuscript 
#This is the v1 version that includes excel processed file

#load packages 
library(tidyverse)
library(ggthemes)
library(scales)
library(ggpubr)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(stats)
library("writexl")
library(PMCMRplus)
library("ggbeeswarm")
library(rstatix)

#-------------------------
#import .xlxs dataset
library(readxl)
adhesioncombined1 <- read_excel("Downloads/adhesion_assay_combo/adhesioncombined1.xlsx")
View(adhesioncombined1)

#---------------------------
#add column that defines adhesion reagent type 

adhesion_combo <- adhesioncombined1 %>% 
  #FIRST ASSIGN DONOR 
  mutate(assay_type = case_when(
    #Assay1
    ScreenID == "3371" ~"CD45+"
    #Assay2
    ,ScreenID == "3391" ~"Manual"
    #Assay3
    ,ScreenID == "3389" ~"BlueWasher"
    #Assay4
    ,ScreenID == "3370" ~"Sunbright"
  )
  )

#View(adhesion_combo)

#Make assay type labels in desired order
adhesion_combo$assay_type <- factor(adhesion_combo$assay_type, levels = c("BlueWasher","Manual","CD45+","Sunbright"))

#-------------------------------------------
#Summarize number of observations for each adhesion reagent type 
adhesion_1 <- adhesion_combo %>%
  group_by(assay_type, nuc_final) %>%
  dplyr::summarise(counts = n()) 

#call number of observations
adhesion_1

#------------------------------------------
#Create violin plot based on dataset  

adhesion_plot <- ggplot(adhesion_combo, aes(x = assay_type, y = nuc_final, fill = assay_type)) +
  geom_violin(show.legend = FALSE)+ 
  #geom_text(data = adhesion_1_label,aes(x=assay_type,y=-20,label=N), size=2) +
  #scale_y_continuous(limits = c(0, 15000)) +
  ylab("Number of recovered cells") + 
  xlab("Assay Type") 

#Run Wilcoxin to compare
compare_means(nuc_final ~ assay_type,  data = adhesion_combo, paired = FALSE)
my_comparisons = list( c("BlueWasher", "Manual"), c("Manual", "CD45+"), c("Sunbright", "Manual"), c("BlueWasher", "CD45+"), c("BlueWasher", "Sunbright"))
stat_compare_means(comparisons = my_comparisons)

#Print violin plot with Wilcoxin results
adhesion_plot + geom_boxplot(width = 0.05, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("BlueWasher" = "BlueWasher", "Manual "= "Manual", "CD45+" = "CD45+", "Sunbright" = "Sunbright")) 


