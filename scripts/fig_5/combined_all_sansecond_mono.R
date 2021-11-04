library(tidyverse)
library(ggthemes)
library(scales)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(stats)
#library(plyr)
library(hexbin)
library(plotly)
library(MASS)
library(patchwork)
library(ggExtra)

#find and read data from the directory
library(readxl)
#setwd("/Users/gallantkl")
#g_all301 <-combined_nosecondary <- read_excel("Downloads/merge_301/combined_nosecondary.xlsx")

g_all301 <- read_excel("data/raw-data_v1/fig_5/combined_nosecondary.xlsx")
#View(g_all301)

#filter CD4 data 
g_Mono <- filter(g_all301, cell_type == "Monocyte") 


g_Mono$treat <- factor(g_Mono$treat, levels = c("Untreated","ETP"))

g_Mono2 <- g_Mono %>%
  mutate(repli = case_when(
    Row %in% c(3,9) ~ "1", 
    Row %in% c(4,10) ~ "2", 
    Row %in% c(5,11) ~ "3"
  ))

View(g_Mono2)

#counts 
m_Mono <- g_Mono2 %>%
  group_by(treat, donor) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

sum_Mono <- g_Mono2 %>%
  group_by(treat, donor) %>%
  summarise(mean = mean(nuc_intyH2AX_mean))


sum_Mono_2 <- g_Mono2 %>%
  group_by(treat, donor) %>%
  summarise(mean = mean(spot_int53BP1_sum))

all_repli <- g_Mono2 %>% 
  group_by(treat,donor, repli) %>% 
  summarise(mean1 = mean(spot_int53BP1_sum), 
            mean2= mean(nuc_intyH2AX_mean))
#all_repli2 <- g_Mono2 %>% 
 # group_by(treat,donor, repli) %>% 
 # summarise(mean = mean(nuc_intyH2AX_mean))


#ANOVA on untreated 

#53BP1 
all_repli_un <- filter (all_repli,treat == "Untreated" )

compare_means(mean1 ~ repli, data = all_repli_un, method = "anova")
#anova
ad_aov_2 <- aov(mean1 ~ as.factor(repli) + as.factor(donor), data = all_repli_un)
summary(ad_aov_2)

#yH2AX 
compare_means(mean2 ~ repli, data = all_repli_un, method = "anova")
#anova
ad_aov_3 <- aov(mean2 ~ as.factor(repli) + as.factor(donor), data = all_repli_un)
summary(ad_aov_3)


#ANOVA on treated

all_repli_ETP <- filter (all_repli,treat == "ETP" )

compare_means(mean1 ~ repli, data = all_repli_ETP, method = "anova")
#anova
ad_aov_4 <- aov(mean1~ as.factor(repli) + as.factor(donor), data = all_repli_ETP)
summary(ad_aov_4)

compare_means(mean2 ~ repli, data = all_repli_ETP, method = "anova")
#anova
ad_aov_5 <- aov(mean2~ as.factor(repli) + as.factor(donor), data = all_repli_ETP)
summary(ad_aov_5)


range_mono_53bp1 <- aggregate(spot_int53BP1_sum ~ donor + treat, g_Mono2, function(g_Mono2) max(g_Mono2[g_Mono2>0], na.rm=TRUE) / min(g_Mono2[g_Mono2>0], na.rm=TRUE) )
range_mono_yH2AX <- aggregate(nuc_intyH2AX_mean ~ donor + treat, g_Mono2, function(g_Mono2) max(g_Mono2[g_Mono2>0], na.rm=TRUE) / min(g_Mono2[g_Mono2>0], na.rm=TRUE) )


range_mono_53bp1$log_base2 = log2(range_mono_53bp1$spot_int53BP1_sum)
range_mono_yH2AX$log_base2 = log2(range_mono_yH2AX$nuc_intyH2AX_mean)

aggregate(log_base2 ~ treat, range_mono_53bp1, mean)
aggregate(log_base2 ~ treat, range_mono_yH2AX, mean)


#TOTAL SCATTERPLOT--------
#total graphs

#make continous variable categorical 
g_Mono2$cat_donor<-cut(g_Mono2$donor, seq(1,11,1), right=FALSE, labels=c(1:10))


#View(donor_untreated)

#created scatterplot with untreated data-> categorical 
g_Mono_all <- ggplot(g_Mono2, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=cat_donor)) + 
  geom_point(alpha = 0.5)+
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug() +
  scale_x_continuous(trans = log2_trans()) + 
  scale_y_continuous(trans = log2_trans())

#run untreated data scatterplot 
g_Mono_all + facet_grid(treat ~ donor) + scale_color_grey(start = 0.3, end = 0.3) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  


#-------------------------------
#contour plot 
g_Mono2$cat_donor<-cut(g_Mono2$donor, seq(1,11,1), right=FALSE, labels=c(1:10))


#View(donor_untreated)

#created scatterplot with untreated data-> categorical 
g_Mono_contour <- ggplot(g_Mono2, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color = treat)) + 
  geom_density2d(bins = 50) +
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  #geom_rug() +
  scale_x_continuous(trans = 'log2') + 
  scale_y_continuous(trans = 'log2')

g_Mono_contour

#run untreated data scatterplot 
g_Mono_contour + facet_wrap(~donor, nrow = 2, ncol = 5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

#_______________________________
#another filled density/contour plot. 

g_Mono_confilled <- ggplot(g_Mono2, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat), show.legend = FALSE) +
  xlab("53BP1 sum spot intensity (a.u.)") + 
  ylab(expression(gamma ~ "H2AX mean nuclear intensity (a.u.)")) + 
  labs(fill="Treatment") + 
  ggtitle("Monocytes") + 
  scale_x_continuous(trans = 'log2', labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific) 

#g_Mono_confilled + facet_wrap(~donor, nrow =2, ncol =5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
 #g_Mono_confilled + facet_wrap(~donor, nrow =2, ncol =5) + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text=element_text(size=15))

 #g_Mono_confilled + facet_wrap(~donor, nrow =2, ncol =5) + theme_classic(base_size = 20) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text=element_text(size=15))
 
 g_Mono_confilled + facet_wrap(~donor, nrow = 2, ncol = 5) + theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=15), panel.grid.major = element_blank(),
                                                                             panel.grid.minor = element_blank(), strip.background = element_blank())

dm2 <- ggMarginal(g_Mono_confilled, type = "density")

dens1 <- ggplot(g_Mono2, aes(x = spot_int53BP1_sum, fill = treat)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")

dens2 <- ggplot(g_Mono2, aes(x = nuc_intyH2AX_mean, fill = treat)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()


dens1 + plot_spacer() + g_Mono_confilled + dens2 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

#________________________________________________________________________________________
kd <- with(g_Mono2, MASS::kde2d(x= spot_int53BP1_sum + 1, y= nuc_intyH2AX_mean + 1))
fig <- plot_ly(x = spot_int53BP1_sum, y = nuc_intyH2AX_mean, z = count) %>% add_surface()

fig


fig <- plot_ly(g_Mono2, x = ~spot_int53BP1_sum, y = ~nuc_intyH2AX_mean) 
fig <- fig %>%
  add_trace(
    type='histogram2dcontour',
    contours = list(
      showlabels = T,
      labelfont = list(
        family = 'Raleway',
        color = 'white'
      )
    ),
    hoverlabel = list(
      bgcolor = 'white',
      bordercolor = 'black',
      font = list(
        family = 'Raleway',
        color = 'black'
      )
    )
  )

fig
#------------------------------------------------------------------------------------

#Quick Scatterplot 
g_Mono_scatter <-ggplot(g_Mono2, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=treat)) + 
  geom_point(alpha = 0.5)+
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug() +
  scale_x_log10() + 
  scale_y_log10() 
#geom_smooth(method= lm, se=FALSE, fullrange=TRUE)  

#etp_CD4 <- filter(m_CD4, treat == "ETP")

#etp_CD4_Reg <- lm(spot_intyH2AX_mean ~ spot_int53BP1_sum, data = etp_CD4) 
#summary(etp_CD4_Reg)

#summary(etp_CD4_Reg)$r.squared



#unt_CD4 <- filter(m_CD4, treat == "Untreated")

#unt_CD4_Reg <- lm(spot_intyH2AX_mean ~ spot_int53BP1_sum, data = unt_CD4) 
#summary(unt_CD4_Reg)

#summary(unt_CD4_Reg)$r.squared



#scatter_CD4 + facet_grid(vars(treat)) + labs(fill="Treatment") 

g_Mono_scatter + labs(fill="treat") + facet_wrap(~ donor, ncol = 2, nrow =5)

#Filter monocyte data to untreated 
donor_untreated <- filter(g_Mono2, treat == "Untreated")
#donor_untreated$c_donor <- cut(donor_untreated$donor, breaks = c(1,2,3,4,5,6,7,8,9,10), 
                             #  labels = c(1,2,3,4,5,6,7,8,9,10))

donor_untreated$cat_donor<-cut(donor_untreated$donor, seq(1,11,1), right=FALSE, labels=c(1:10))

#donor_untreated$c_donor <-revalue(donor_untreated$donor, c( "1" = "A", "2" = "B", "3" = "C", "4" = "D", "5" = "E", "6" = "F", "7" = "G", "8" = "H", "9" = "I", "10" ="J"))

View(donor_untreated)

#created scatterplot with untreated data-> categorical 
g_Mono_untreat <- ggplot(donor_untreated, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=cat_donor)) + 
  geom_point(alpha = 0.5)+
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug() +
  scale_x_continuous(trans = log2_trans()) + 
  scale_y_continuous(trans = log2_trans())  
 
  

#create scatterplot non-categorical 
g_Mono_untreat <- ggplot(donor_untreated, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=donor)) + 
  geom_point(alpha = 0.5)+
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug() +
  scale_x_log10() + 
  scale_y_log10() + 
  
  
g_Mono_untreat

#run untreated data scatterplot 
g_Mono_untreat + facet_wrap(~ donor, ncol = 5, nrow =2) + scale_color_grey(start = 0.1, end = 0.1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

#Filter monocyte data to untreated 
donor_ETP <- filter(g_Mono2, treat == "ETP")
#donor_untreated$c_donor <- cut(donor_untreated$donor, breaks = c(1,2,3,4,5,6,7,8,9,10), 
#  labels = c(1,2,3,4,5,6,7,8,9,10))

donor_ETP$cat_donor<-cut(donor_ETP$donor, seq(1,11,1), right=FALSE, labels=c(1:10))

#donor_untreated$c_donor <-revalue(donor_untreated$donor, c( "1" = "A", "2" = "B", "3" = "C", "4" = "D", "5" = "E", "6" = "F", "7" = "G", "8" = "H", "9" = "I", "10" ="J"))

View(donor_ETP)

#created scatterplot with treated data-> categorical 
g_Mono_ETP <- ggplot(donor_ETP, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=cat_donor)) + 
  geom_point(alpha = 0.5)+
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug() +
  scale_x_continuous(trans = log2_trans()) + 
  scale_y_continuous(trans = log2_trans())  


#create scatterplot non-categorical 
g_Mono_ETP <- ggplot(donor_ETP, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=density)) + 
  geom_point(alpha = 0.5)+
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug() +
  scale_x_log10() + 
  scale_y_log10() 

#run treated data scatterplot 
g_Mono_ETP + facet_wrap(~ donor, ncol = 5, nrow =2) + scale_color_grey() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#spots 53BP1
g_53BP1 <- ggplot(g_Mono2, aes(x= treat, y= spots_53BP1 +1, fill= treat)) +
  geom_violin() + 
  scale_y_continuous(trans = "log10", limits = c(1e0, 1e+1)) +
  stat_summary(fun.data = n_fun, geom = "text") +
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("CD4") +
  labs(fill="Treatment") +
  theme(axis.line = element_line(colour="black", size = 3), panel.border = element_rect(colour = "black", fill=NA, size =2))

#print violin plot 
g_53BP1  + geom_boxplot(width=0.1, notch = TRUE, boxlwd = 2) + scale_fill_brewer(palette = "BuGn") + theme_bw() + facet_wrap(~donor, ncol = 2, nrow =3)




gbar_Mono <- ggplot(g_Mono2, aes(x = treat, y = spots_53BP1)) +
  geom_bar() + 
  ggtitle("CD4") + 
  labs(x = "Treatment", y = "53BP1 Spots") 

gbar_Mono + facet_wrap(~donor)

gbar_Mono


#density plot attempt 
dense_mono <- ggplot(donor_untreated, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) +
  geom_bin2d(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

dense_mono


dense_mono + facet_wrap(~ donor, ncol = 2, nrow =5)



#densityplot trial 2

g_Mono_dense <- ggplot(g_Mono2, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=density)) + 
  geom_point(alpha = 0.5)+
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") 
  
g_Mono_dense


g_Mono_dense + facet_wrap(~ donor, ncol = 2, nrow =5)


g_Mono_test <- ggplot(donor_untreated, aes(x = spot_int53BP1_sum, y = donor)) +
  geom_density_ridges(aes(fill = donor)) 
 

#-----------------------------------------------------------------
#density plot 
g_Mono_hex <- ggplot(g_Mono2, mapping = aes(x= spot_int53BP1_sum , y= nuc_intyH2AX_mean)) + 
  geom_hex(bins = 50) +
  scale_fill_viridis_c() + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  scale_x_continuous(trans = 'log2') + 
  scale_y_continuous(trans = 'log2')

g_Mono_hex

#run geom hexbin
g_Mono_hex + facet_grid(treat ~ donor) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
