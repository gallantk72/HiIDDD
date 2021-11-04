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
#setwd("/Users/gallantkl")
#import dataset
library(readxl)
#adhesioncombined1 <- read_excel("Downloads/adhesion_assay_combo/adhesioncombined1.xlsx")


adhesioncombined1 <- read_excel("data/raw-data_v1/fig_1/fig_1b/adhesioncombined1.xlsx")
View(adhesioncombined1)

#-------------------------------------------------------------------------
#prepare dataset for different cell types 

#Add column to define assay type 

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

View(adhesion_combo)

adhesion_combo$assay_type <- factor(adhesion_combo$assay_type, levels = c("BlueWasher","Manual","CD45+","Sunbright"))


#------------------------------------

#Total counts for each well 
adhesion_1 <- adhesion_combo %>%
  group_by(assay_type, nuc_final) %>%
  dplyr::summarise(counts = n()) 

adhesion_2 <- adhesion_combo %>%
  group_by(assay_type) %>%
  summarise(mean = mean(nuc_final)

#Per well data 
adhesion_1_label <- adhesion_combo %>%
  group_by(assay_type, nuc_final) %>%
  summarise(counts = n()) %>%
  summarise(N=paste0("n =", n()))


aggregate(nuc_final ~ assay_type, adhesion_combo,min)
aggregate(WellName ~ assay_type, adhesion_combo, counts)

adhesion_plot <- ggplot(adhesion_combo, aes(x = assay_type, y = nuc_final, fill = assay_type)) +
  #geom_dotplot(aes(fill = assay_type), binaxis = "y", binwidth = 0.1, stackdir= "center" ) + 
  geom_violin(show.legend = FALSE)+ 
  #geom_jitter()+ 
  #geom_text(data = adhesion_1_label,aes(x=assay_type,y=-20,label=N), size=2) +
  #scale_y_continuous(limits = c(0, 15000)) +
  ylab("Number of recovered cells") + 
  xlab("Assay Type") 

compare_means(nuc_final ~ assay_type,  data = adhesion_combo, paired = FALSE)
my_comparisons = list( c("BlueWasher", "Manual"), c("Manual", "CD45+"), c("CD45+", "Sunbright"), c("Sunbright", "Manual"), c("BlueWasher", "CD45+"), c("BlueWasher", "Sunbright"))
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_final ~ assay_type, data = adhesion_combo)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(adhesion_combo$nuc_final,   adhesion_combo$assay_type, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

aggregate(nuc_final ~ assay_type, adhesion_combo,mean)
aggregate(nuc_final ~ assay_type, adhesion_combo,sd)
#scale_x_discrete(labels=c("0.5" = "Dose 0.5", "1" = "Dose 1",
#"2" = "Dose 2"))
# look at the comparisons and p-values
head(tidy_ad_pairwise)
#print violin plot with stats 
adhesion_plot + geom_boxplot(width = 0.05, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#+ scale_x_discrete(labels = c("BlueWasher"= "BlueWasher n=57", "Manual"= "Manual n=64", "60000"="6", "CD45+"= "CD45+ n=47", "Sunbright" = "Sunbright n=46")) 

adhesion_plot + geom_boxplot(width = 0.05, show.legend = FALSE) +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) 
