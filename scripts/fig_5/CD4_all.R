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

library(rstatix)

#import dataset
#setwd("/Users/gallantkl")
library(readxl)
#CD4_combined <- read_excel("Documents/CD4_training/CD4_combined.xlsx")
CD4_combined <- read_excel("data/raw-data_v1/fig_5/CD4_combined.xlsx")
View(CD4_combined)

#CD4_combined %>% 
  #mutate(donor = case_when(
    #Donor1
    #starts_with(WellName %in% c('C','D','E') & ends_with(WellName %in% c('3', '5'))) ~"1"
    #Donor2
    #,starts_with(WellName %in% c(C,D,E) & ends_with(WellName %in% c('7', '9'))) ~"2"
    #Donor 3 
    #,starts_with(WellName %in% c(C,D,E) & ends_with(WellName %in% c('11', '13'))) ~"3"
    #Donor 4
    #,starts_with(WellName %in% c(C,D,E) & ends_with(WellName %in% c('15', '17'))) ~"4"
    #Donor 5
    #,starts_with(WellName %in% c(C,D,E) & ends_with(WellName %in% c('19', '21'))) ~"5"
    #Donor 6
   # ,starts_with(WellName %in% c(K,L,M) & ends_with(WellName %in% c('3', '5'))) ~"6"
    #Donor 7
    #,starts_with(WellName %in% c(K,L,M) & ends_with(WellName %in% c('7', '9'))) ~"7"
    #Donor 8 
    #,starts_with(WellName %in% c(K,L,M) & ends_with(WellName %in% c('11', '13'))) ~"8"
    #Donor 9 
    #,starts_with(WellName %in% c(K,L,M) & ends_with(WellName %in% c('15', '17'))) ~ "9"
    #Donor 10
    #,starts_with(WellName %in% c(K,L,M) & ends_with(WellName %in% c('19', '21'))) ~"10"
#  ))


#CD4_combined %>% 
  #mutate(donor = case_when(
    #Donor1
   # starts_with(WellName %in% c(C,D,E)) & ends_with(WellName %in% c("3", "5"))) ~"1"
 # )

#-------------------------------------------------
#Prepare dataset 

CD4_combined2 <- CD4_combined %>% 
  #FIRST ASSIGN DONOR 
  mutate(donor = case_when(
    #Donor1
    WellName %in% c("C3","C5","D3", "D5", "E3", "E5") ~"1"
    #Donor2
    ,WellName %in% c("C7","C9","D7", "D9", "E7", "E9") ~"2"
    #Donor3
    ,WellName %in% c("C11","C13","D11", "D13", "E11", "E13") ~"3"
    #Donor5
    ,WellName %in% c("C15","C17","D15", "D17", "E15", "E17") ~"5"
    #Donor4
    ,WellName %in% c("C19","C21","D19", "D21", "E19", "E21") ~"4"
    #Donor6
    ,WellName %in% c("K3","K5","L3", "L5", "M3", "M5") ~"6"
    #Donor7
    ,WellName %in% c("K7","K9","L7", "L9", "M7", "M9") ~"7"
    #Donor8
    ,WellName %in% c("K11","K13","L11", "L13", "M11", "M13") ~"8"
    #Donor9
    ,WellName %in% c("K15","K17","L15", "L17", "M15", "M17") ~"9"
    #Donor10
    ,WellName %in% c("K19","K21","L19", "L21", "M19", "M21") ~"10")
    #ASSIGN TREATMENT COLUMN 
    ,treat = case_when(
      #ETP
      Column %in% c("5","9","13", "17", "21") ~"ETP"
      #Untreated
      ,Column %in% c("3","7","11", "15", "19") ~"Untreated"),
    #ASSIGN replicates
    repli = case_when(
      Row %in% c("3","11") ~ "1", 
      Row %in% c("4","12") ~ "2", 
      Row %in% c("5","13") ~ "3"
    )
    )


#assign treatment
#CD4_combined2 <- CD4_combined %>% 
  #mutate(treat = case_when(
    #Donor1
    #Column %in% c("5","9","13", "17", "21") ~"ETP"
    #,Column %in% c("3","7","11", "15", "19") ~"Untreated"
    #))

#View new dataframe
View(CD4_combined2)

#Export dataframe
write_xlsx(CD4_combined2,"C:\\Users\\gallantkl\\Documents\\CD4_training\\CD4_combined2.xlsx")

CD4_combined2$treat <- factor(CD4_combined2$treat, levels = c("Untreated","ETP"))
#----------------------------------------
#Summarize data by values 

sum_CD4 <- CD4_combined2 %>%
  group_by(treat, donor) %>%
  summarise(mean = mean(nuc_intyH2AX_mean))

View(sum_CD4)


sum_CD4_2 <- CD4_combined2 %>%
  group_by(treat, donor) %>%
  summarise(mean = mean(spot_int53BP1_sum))

View(Sum_CD4_2)

#counts 
#counts 
m_CD4 <- CD4_combined2 %>%
  group_by(treat, donor) %>%
  dplyr::summarise(counts = n())


%>%
  mutate(perc=counts/sum(counts))

#Summarize data by donor and replicates 
donor_repli <-CD4_combined2 %>%
  group_by(treat,donor,repli) %>% 
  summarise(mean = mean(spot_int53BP1_sum))

View(donor_repli)

donor_repli_2 <- CD4_combined2 %>%
  group_by(treat,donor) %>% 
  get_summary_stats(spot_int53BP1_sum, type="mean_sd") %>%
  filter(treat == "Untreated")

View(donor_repli_2)

repli_CD4 <- ggplot(donor_repli_2, aes(x = donor, y = mean)) +
  geom_dotplot(aes(fill = repli), binaxis = "y", binwidth = 0.1, stackdir= "center" ) + 
  stat_summary(fun.y = mean(), fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.1) +
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) + 
  #scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  xlab("Treatment") + 
  #labs(fill="Treatment") + 
  ggtitle("CD4") + 
  theme_bw()

compare_means(mean ~ repli + donor, data = donor_repli_2, method = "anova") 
#------------------------------------------------------------

CD4_combined2$donor <- factor(CD4_combined2$donor, levels = c("1","2", "3", "4", "5", "6", "7", "8", "9", "10"))
#Create contour plot 
CD4_combined2_confilled <- ggplot(CD4_combined2, aes(x= spot_int53BP1_sum , y= nuc_intyH2AX_mean)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat),show.legend = FALSE ) +
  xlab("53BP1 sum spot intensity (a.u.)") + 
  ylab(expression(gamma ~ "H2AX mean nuclear intensity (a.u.)")) + 
  labs(fill="Treatment") + 
  ggtitle("CD4+T cells") + 
  scale_x_continuous(trans = 'log2', labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific) 

#CD4_combined2_confilled + facet_wrap(~donor, nrow =2, ncol =5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
#CD4_combined2_confilled + facet_wrap(~donor, nrow = 2, ncol = 5) + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text=element_text(size=15))  
#CD4_combined2_confilled + facet_wrap(~donor, nrow = 2, ncol = 5) + theme_classic(base_size = 20) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text=element_text(size=15), )  

CD4_combined2_confilled + facet_wrap(~donor, nrow = 2, ncol = 5) + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=15), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), strip.background = element_blank())
#CD4_combined2_confilled + facet_wrap(~donor)
#facet_wrap(~donor, nrow = 2, ncol = 5)

cd4_untreat <- filter(CD4_combined2, treat == "Untreated")
cd4_treat <- filter(CD4_combined2, treat == "ETP")
#find range of variable x and y in the data frame
sapply(CD4_combined2[c('spot_int53BP1_sum','nuc_intyH2AX_mean')], function(CD4_combined2) max(CD4_combined2, na.rm=TRUE) - min(CD4_combined2, na.rm=TRUE))
sapply(cd4_untreat[c('spot_int53BP1_sum','nuc_intyH2AX_mean')], function(cd4_untreat) max(cd4_untreat, na.rm=TRUE) / min(cd4_untreat, na.rm=TRUE))
sapply(cd4_untreat[c('spot_int53BP1_sum','nuc_intyH2AX_mean')], function(cd4_untreat) min(cd4_untreat, na.rm=TRUE))
sapply(cd4_untreat[c('spot_int53BP1_sum','nuc_intyH2AX_mean')], function(cd4_untreat) max(cd4_untreat, na.rm=TRUE))

sapply(cd4_treat[c('spot_int53BP1_sum','nuc_intyH2AX_mean')], function(cd4_treat) max(cd4_treat, na.rm=TRUE) / min(cd4_treat, na.rm=TRUE))
sapply(cd4_treat[c('spot_int53BP1_sum','nuc_intyH2AX_mean')], function(cd4_treat) min(cd4_treat, na.rm=TRUE))
sapply(cd4_treat[c('spot_int53BP1_sum','nuc_intyH2AX_mean')], function(cd4_treat) max(cd4_treat, na.rm=TRUE))
#find range of all variables in the data frame
sapply(df, function(df) max(df, na.rm=TRUE) - min(df, na.rm=TRUE))

#-----------------------------------------------------------------
#ANOVA 

all_repli_CD4 <- CD4_combined2 %>% 
  group_by(treat,donor, repli) %>% 
  summarise(mean1 = mean(spot_int53BP1_sum), 
            mean2= mean(nuc_intyH2AX_mean))
#all_repli_CD42 <- g_Mono2 %>% 
# group_by(treat,donor, repli) %>% 
# summarise(mean = mean(nuc_intyH2AX_mean))


cd4_repli <- CD4_combined2 %>% 
  group_by(treat, donor)

cd4_repli_min <- CD4_combined2 %>%
  group_by(donor) %>%
  get_summary_stats(spot_int53BP1_sum, type = c("min"))

cd4_repli_max <- CD4_combined2 %>%
  group_by(donor) %>%
  get_summary_stats(spot_int53BP1_sum, type = c("max"))

cd4_repli_min2 <- CD4_combined2 %>%
  group_by(donor) %>%
  get_summary_stats(nuc_intyH2AX_mean, type = c("min"))

cd4_repli_max2 <- CD4_combined2 %>%
  group_by(donor) %>%
  get_summary_stats(nuc_intyH2AX_mean, type = c("max"))

cd4_repli_max - cd4_repli_min   

cd4_repli <- CD4_combined2 %>%
  get_summary_stats(spot_int53BP1_sum, type = c("min","max"))
  
 sapply(cd4_repli[c('spot_int53BP1_sum', 'nuc_intyH2AX_mean')], function(cd4_repli) max(cd4_repli[cd4_repli>0], na.rm=TRUE) / min(cd4_repli[cd4_repli>0], na.rm=TRUE))#ANOVA on untreated 

range_cd4_53bp1 <- aggregate(spot_int53BP1_sum ~ donor + treat, CD4_combined2, function(cd4_repli) max(cd4_repli[cd4_repli>0], na.rm=TRUE) / min(cd4_repli[cd4_repli>0], na.rm=TRUE) )
range_cd4_yH2AX <- aggregate(nuc_intyH2AX_mean ~ donor + treat, CD4_combined2, function(cd4_repli) max(cd4_repli[cd4_repli>0], na.rm=TRUE) / min(cd4_repli[cd4_repli>0], na.rm=TRUE) )


range_cd4_53bp1$log_base2 = log2(range_cd4_53bp1$spot_int53BP1_sum)
range_cd4_yH2AX$log_base2 = log2(range_cd4_yH2AX$nuc_intyH2AX_mean)

aggregate(log_base2 ~ treat, range_cd4_53bp1, mean)
aggregate(log_base2 ~ treat, range_cd4_yH2AX, mean)

#53BP1 
all_repli_CD4_un <- filter (all_repli_CD4,treat == "Untreated" )

compare_means(mean1 ~ repli, data = all_repli_CD4_un, method = "anova")
#anova
ad_aov_2 <- aov(mean1 ~ as.factor(repli) + as.factor(donor), data = all_repli_CD4_un)
summary(ad_aov_2)

TukeyHSD(ad_aov_2)

#yH2AX 
compare_means(mean2 ~ repli, data = all_repli_CD4_un, method = "anova")
#anova
ad_aov_3 <- aov(mean2 ~ as.factor(repli) + as.factor(donor), data = all_repli_CD4_un)
summary(ad_aov_3)

TukeyHSD(ad_aov_3)

#ANOVA on treated

all_repli_CD4_ETP <- filter (all_repli_CD4,treat == "ETP" )

#53BP1
compare_means(mean1 ~ repli, data = all_repli_CD4_ETP, method = "anova")
#anova
ad_aov_4 <- aov(mean1~ as.factor(repli) + as.factor(donor), data = all_repli_CD4_ETP)
summary(ad_aov_4)

#yH2AX
compare_means(mean2 ~ repli, data = all_repli_CD4_ETP, method = "anova")
#anova
ad_aov_5 <- aov(mean2~ as.factor(repli) + as.factor(donor), data = all_repli_CD4_ETP)
summary(ad_aov_5)
