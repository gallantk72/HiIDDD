library(tidyverse)
library(ggthemes)
library(scales)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(stats)
library(plyr)

#find and view dataset from directory 


library(readxl)
#setwd("/Users/gallantkl")
#library(readxl)
#CD8_E13 <- read_excel("Downloads/merge_302/combined_302_excel.xls")
#View(combined_302_excel)

library(readxl)
#CD8_E13 <- read_excel("Downloads/merge_302/CD8_all_sansE13.xlsx")
CD8_E13 <- read_excel("data/raw-data_v1/fig_5/CD8_all_sansD13.xlsx")
View(CD8_E13)


#Re-order groups for data  
CD8_E13$treat <- factor(CD8_E13$treat, levels = c("Untreated","ETP"))

#counts 
m_CD8 <- CD8_E13 %>%
  group_by(treat, donor) %>%
  dplyr::summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))


sum_CD8 <- CD8_E13 %>%
  group_by(treat, donor) %>%
  summarise(mean = mean(nuc_intyH2AX_mean))


sum_CD8_2 <- CD8_E13 %>%
  group_by(treat, donor) %>%
  summarise(mean = mean(spot_int53BP1_sum))

THP1_1perc <- THP1_1 %>%
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))
 

aggregate(spot_int53BP1_sum ~ treat, j_fil_Jur,mean)
aggregate(spot_int53BP1_sum ~ treat, j_fil_Jur,sd)
aggregate(nuc_intyH2AX_mean ~ treat, j_fil_Jur,mean)
aggregate(nuc_intyH2AX_mean ~ treat, j_fil_Jur,sd)

#____________________________________________
#total graphs

#make continous variable categorical 
CD8_E13$cat_donor<-cut(CD8_E13$donor, seq(1,11,1), right=FALSE, labels=c(1:10))


#View(donor_untreated)

#created scatterplot with untreated data-> categorical 
CD8_E13_all <- ggplot(CD8_E13, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=cat_donor)) + 
  geom_point(alpha = 0.5)+
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug() +
  scale_x_continuous(trans = log2_trans()) + 
  scale_y_continuous(trans = log2_trans())

#run untreated data scatterplot 
CD8_E13_all + facet_grid(treat ~ donor) + scale_color_grey(start = 0.3, end = 0.3) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text=element_text(size=15))



#__________________________________________________________________________________________________
range_cd8_53bp1 <- aggregate(spot_int53BP1_sum ~ donor + treat, CD8_E13, function(CD8_E13) max(CD8_E13[CD8_E13>0], na.rm=TRUE) / min(CD8_E13[CD8_E13>0], na.rm=TRUE) )
range_cd8_yH2AX <- aggregate(nuc_intyH2AX_mean ~ donor + treat, CD8_E13, function(CD8_E13) max(CD8_E13[CD8_E13>0], na.rm=TRUE) / min(CD8_E13[CD8_E13>0], na.rm=TRUE) )


range_cd8_53bp1$log_base2 = log2(range_cd8_53bp1$spot_int53BP1_sum)
range_cd8_yH2AX$log_base2 = log2(range_cd8_yH2AX$nuc_intyH2AX_mean)

aggregate(log_base2 ~ treat, range_cd8_53bp1, mean)
aggregate(log_base2 ~ treat, range_cd8_yH2AX, mean)
#geom_hexbin 

CD8_E13_hex <- ggplot(CD8_E13, mapping = aes(x= spot_int53BP1_sum , y= nuc_intyH2AX_mean)) + 
  geom_hex(bins = 60) +
  scale_fill_viridis_c() + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  scale_x_continuous(trans = 'log2') + 
  scale_y_continuous(trans = 'log2')

CD8_E13_hex

#run geom hexbin
CD8_E13_hex + facet_grid(treat ~ donor) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#--------
#geom_hexbar with both (facet_warap) 

#geom_hexbin 

CD8_E13_hex <- ggplot(CD8_E13, mapping = aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) + 
  geom_hex(bins = 20) +
  scale_fill_viridis_c() + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  scale_x_continuous(trans = log2_trans()) + 
  scale_y_continuous(trans = log2_trans())

CD8_E13_hex

#run untreated data scatterplot 
CD8_E13_hex + facet_wrap(~ donor, ncol = 5, nrow =2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  


#-----------
#contour plot 
CD8_E13$cat_donor<-cut(CD8_E13$donor, seq(1,11,1), right=FALSE, labels=c(1:10))


#View(donor_untreated)

#created scatterplot with untreated data-> categorical 
CD8_E13_contour <- ggplot(CD8_E13, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color = treat)) + 
  geom_density2d(bins = 50) +
  xlab("53BP1 Sum Spot Intensity (A.U.)") + 
  ylab("yH2AX Mean Nuclear Intensity (A.U.)") + 
  labs(fill="Treatment") + 
  #geom_rug() +
  scale_x_continuous(trans = 'log2') + 
  scale_y_continuous(trans = 'log2')

CD8_E13_contour

#run untreated data scatterplot 
#CD8_E13_contour + facet_wrap(~donor, nrow = 2, ncol = 5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
CD8_E13_contour + facet_wrap(~donor, nrow = 2, ncol = 5) + theme_bw(base_size = 20)  
#_______________________________
#another filled density/contour plot. 

CD8_E13_confilled <- ggplot(CD8_E13, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat),show.legend = FALSE ) +
  xlab("53BP1 sum spot intensity (a.u.)") + 
  ylab(expression(gamma ~ "H2AX mean nuclear intensity (a.u.)")) + 
  labs(fill="Treatment", titles = expression("CD8+T cells")) + 
  #ggtitle("CD8+T cells") + 
  scale_x_continuous(trans = 'log2', labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific)  
 

#CD8_E13_confilled + facet_wrap(~donor, nrow =2, ncol =5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
#CD8_E13_confilled + facet_wrap(~donor, nrow =2, ncol =5) + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text=element_text(size=15))  
#CD8_E13_confilled + facet_wrap(~donor, nrow =2, ncol =5) + theme_classic(base_size = 20) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text=element_text(size=15))  

CD8_E13_confilled + facet_wrap(~donor, nrow = 2, ncol = 5) + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=15), panel.grid.major = element_blank(),
                                                                                                    panel.grid.minor = element_blank(), strip.background = element_blank())
#--------
#another filled density plot for no of 53BP1 
CD8_E13_confilled2 <- ggplot(CD8_E13, aes(x= no_pos53BP1+1, y= nuc_intyH2AX_mean+1)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat)) +
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") +   
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10')

CD8_E13_confilled2 + facet_wrap(~donor, nrow =3, ncol =4) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

#---------------------------

CD8_E13$cat_donor<-cut(CD8_E13$donor, seq(1,11,1), right=FALSE, labels=c(1:10))


#View(donor_untreated)

#created scatterplot with untreated data-> categorical 
CD8_E13_contour2 <- ggplot(CD8_E13, aes(x= no_pos53BP1+1, y= nuc_intyH2AX_mean+1, color = treat)) + 
  geom_density2d(bins = 20) +
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  #geom_rug() +
  scale_x_continuous(trans = 'log2') + 
  scale_y_continuous(trans = 'log2')

CD8_E13_contour2

#run untreated data scatterplot 
CD8_E13_contour2 + facet_wrap(~donor, nrow = 2, ncol = 5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

#----------------------------------------------------------------------------------------


#Filter CD8 data to untreated 
donor_untreated <- filter(CD8_E13, treat == "Untreated")

#make continous variable categorical 
donor_untreated$cat_donor<-cut(donor_untreated$donor, seq(1,6,1), right=FALSE, labels=c(1:5))


#View(donor_untreated)

#created scatterplot with untreated data-> categorical 
CD8_E13_untreat <- ggplot(donor_untreated, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=cat_donor)) + 
  geom_point(alpha = 0.5)+
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug() +
  scale_x_continuous(trans = log2_trans()) + 
  scale_y_continuous(trans = log2_trans())



#run untreated data scatterplot 
CD8_E13_untreat + facet_wrap(~ donor, ncol = 10, nrow =1) + scale_color_brewer(palette ="Paired") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

#_________________________________________________________________________________________________________





#Filter CD8 data to treated 
donor_ETP <- filter(CD8_E13, treat == "ETP")

#make continous variable categorical 
donor_ETP$cat_donor<-cut(donor_ETP$donor, seq(1,11,1), right=FALSE, labels=c(1:10))

#make scatterplot
CD8_E13_ETP <- ggplot(donor_ETP, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=cat_donor)) + 
  geom_point(alpha = 0.5)+
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  labs(fill="Treatment") + 
  geom_rug() +
  scale_x_continuous(trans = log2_trans()) + 
  scale_y_continuous(trans = log2_trans())

#run scatterplot
CD8_E13_ETP + facet_wrap(~ donor, ncol = 10, nrow =2) + scale_color_brewer(palette ="Paired") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#-------------------------------------------------------------------
  
  CD8_E13d_ETP <- ggplot(donor_ETP, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) +
  geom_bin2d(bins = 70) + 
  scale_fill_continuous(type = "viridis") + 
  theme_bw()

#-------------------------------------------------
#ANOVA 
#Prepare data

CD8_E13_2 <- CD8_E13 %>%
  mutate(repli = case_when(
    Row %in% c(3,11) ~ "1", 
    Row %in% c(4,12) ~ "2", 
    Row %in% c(5,13) ~ "3"
  ))

View(CD8_E13_2)

#ANOVA 

all_repli_CD8 <- CD8_E13_2 %>% 
  group_by(treat,donor,repli) %>% 
  summarise(mean1 = mean(spot_int53BP1_sum), 
            mean2= mean(nuc_intyH2AX_mean))

View(all_repli_CD8)
#all_repli_CD82 <- g_Mono2 %>% 
# group_by(treat,donor, repli) %>% 
# summarise(mean = mean(nuc_intyH2AX_mean))


#ANOVA on untreated 

#53BP1 
all_repli_CD8_un <- filter (all_repli_CD8,treat == "Untreated" )

compare_means(mean1 ~ repli, data = all_repli_CD8_un, method = "anova")
#anova
ad_aov_2 <- aov(mean1 ~ as.factor(repli) + as.factor(donor), data = all_repli_CD8_un)
summary(ad_aov_2)

TukeyHSD(ad_aov_2)

#yH2AX 
compare_means(mean2 ~ repli, data = all_repli_CD8_un, method = "anova")
#anova
ad_aov_3 <- aov(mean2 ~ as.factor(repli) + as.factor(donor), data = all_repli_CD8_un)
summary(ad_aov_3)

TukeyHSD(ad_aov_3)

#ANOVA on treated

all_repli_CD8_ETP <- filter (all_repli_CD8,treat == "ETP" )

#53BP1
compare_means(mean1 ~ repli, data = all_repli_CD8_ETP, method = "anova")
#anova
ad_aov_4 <- aov(mean1~ as.factor(repli) + as.factor(donor), data = all_repli_CD8_ETP)
summary(ad_aov_4)

#yH2AX
compare_means(mean2 ~ repli, data = all_repli_CD8_ETP, method = "anova")
#anova
ad_aov_5 <- aov(mean2~ as.factor(repli) + as.factor(donor), data = all_repli_CD8_ETP)
summary(ad_aov_5)
