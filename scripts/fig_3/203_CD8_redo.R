library(tidyverse)
library(ggthemes)
library(scales)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(stats)

#find and read data from the directory
library(readxl)
#X190810_203 <- read_excel("Desktop/exp203_CD4_CD8_B_Monocytes_KG/data_203/190810-203-ALL-revised.xlsx")
View(X190810_203)

X190810_203 <- read_excel("data/raw-data_v1/fig_3/190810-203-ALL-revised.xlsx")

#filter CD8 data 
m_CD8 <- filter(X190810_203, cell_type == "CD8") 

#Order two treatment groups 
m_CD8$treat <- factor(m_CD8$treat, levels = c("Untreated","ETP"))


#-----------------------------------------

#Define median and standard deviation for yH2AX nuclear intensity 
aggregate(nuc_intyH2AX_mean ~ treat, m_CD8,mean)
aggregate(nuc_intyH2AX_mean ~ treat, m_CD8,sd)

#filter data for only yH2AX high spots after changing column to threshold 
High_spots_yH2AX_CD8 <- filter(m_CD8,treat_yH2AX == "High")

#Show counts for data with "high" spots 

#Define median and standard deviation for 52BP1 spot intenstiy 
aggregate(spot_int53BP1_sum ~ treat, m_CD8,mean)
aggregate(spot_int53BP1_sum ~ treat, m_CD8,sd)

#filter data for only 53BP1 high spots
treat_lev_CD8 <- filter(m_CD8, treat_lev =="High")

#------------------------------------
#Create bargraphs for CD8 

#Filter each treatment group and treatment elvel for CD8 
#percentage filter 
m_CD8_perc <- m_CD8 %>%
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#filter out percentage of high cells
perc_pos5 <-filter(m_CD8_perc, treat_lev == "High")

#plot barplot based on percentages
bar_CD8 <- ggplot(m_CD8_perc, aes(x = treat, y = perc*100)) +
  geom_bar(
    aes(fill = factor(treat_lev, levels=c("Low","High"))),
    stat = "identity", position = position_stack()
  ) + 
  ggtitle("CD8+T cells") + 
  labs(x = "Treatment", y = "Total Cells (%)", fill = "53BP1 Spot Status") 


bar_CD8 + scale_fill_brewer(palette = "BuGn") + theme_classic()

#--------------------------------------
#create violin plot for 53BP1 

#assign count number for each treatment group
labels_CD8 <-m_CD8 %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))

#this is a sample line of code for how to add # of observations to plot 
#geom_text(data = labels,aes(x=treat,y=5e+5,label=N), size=6)  
  
#create violin plot for 53BP1 
v_53BP1_CD8 <- ggplot(m_CD8, aes(x=treat, y = spot_int53BP1_sum +1, fill = treat)) + 
  geom_violin(size=0.75, show.legend = FALSE)+ 
  scale_y_continuous(trans= "log10", limits = c(1e+0, 1e+8)) + 
  geom_text(data = labels_CD8,aes(x=treat,y=5e+6,label=N), size=6) + 
  ylab("53BP1 spot intensity sum (a.u) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("CD8+T cells") +
  annotate("text", x=1.42, y=50000, label = "61.8%", size = 6) +  #untreated and high
  annotate("text", x=1.42, y=10, label = "38.2%", size = 6) + #untreated and low 
  annotate("text", x=2.43, y=50000, label = "71.2%", size = 6) + #ETP and high 
  annotate("text", x=2.43, y=10, label = "28.8%", size = 6)  #EtP and low 

compare_means(spot_int53BP1_sum ~ treat,  data = treat_lev_CD8, paired = FALSE)
my_comparisons = list( c("ETP", "Untreated") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(spot_int53BP1_sum ~ treat, data = treat_lev_CD8)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(treat_lev_CD8$spot_int53BP1_sum,    treat_lev_CD8$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)


#print violin plot 
v_53BP1_CD8 + stat_compare_means(comparisons = my_comparisons, label = "p.signif") + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + theme(axis.text = element_text(size=15))

fig3b_CD8 <- v_53BP1_CD8 +  scale_fill_brewer(palette = "BuGn") + 
  stat_compare_means(comparisons = my_comparisons, data = treat_lev_CD8, label = "p.signif") + 
  theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) +
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank())

#---------------------------------------------------
#create violin plot for pos-yH2AX

#Create number of observations for pos-yH2AX 
labels_CD8_2 <-High_spots_yH2AX_CD8 %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))

#create violin plot for yH2AX 
v_yH2AX_CD8 <- ggplot(High_spots_yH2AX_CD8, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(size = 0.75, show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels_CD8_2,aes(x=treat,y=5e+5,label=N), size=6) + 
  ylab("Pos-yH2AX Nuclear Intensity Mean (A.U) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("CD8+ T lymphocytes")  

#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = High_spots_yH2AX_CD8, paired = FALSE)
my_comparisons = list( c("ETP", "Untreated") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = High_spots_yH2AX_CD8)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(High_spots_yH2AX_CD8$nuc_intyH2AX_mean,    High_spots_yH2AX_CD8$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")

#print violin plot with boxplot and significance values 
v_yH2AX_CD8 + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) 



#------------------------------------------------------------------

#labels 
#Create number of observations for pos-yH2AX 
labels_CD8 <-m_CD8 %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))


v_yH2AX <- ggplot(m_CD8, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels_CD8,aes(x=treat,y=2e+5,label=N), size=6) + 
  ylab(expression(gamma ~ "H2AX nuclear intensity mean (a.u) + 1")) + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("CD8+T cells")  

#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = m_CD8, paired = FALSE)
my_comparisons = list( c("ETP", "Untreated") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = m_CD8)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(m_CD8$nuc_intyH2AX_mean,    m_CD8$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")

#print violin plot with boxplot and significance values 
v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) 

fig3c_CD8 <- v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + 
  scale_fill_brewer(palette = "BuGn") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) +
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank())
#-------------------------------------
#CHI SQUARE FOR yH2AX
#percentages 
j_CD8_perc <- m_CD8 %>%
  group_by(treat,treat_yH2AX ) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#CHI SQUARE ANALYSIS 
obsfreq <- c(2692,7)
nullprobs <- c(0.5,0.5)
#chisq.test(obsfreq,p=nullprobs)

obsfreq2 <- c(1423,223)
nullprobs <- c(0.5,0.5)

rows = 2 
matriz_CD8 = matrix(c(obsfreq,obsfreq2), 
                    nrow=rows,
                    byrow=TRUE)

rownames(matriz_CD8) = c("DMSO","ETP")
colnames(matriz_CD8) = c("Low", "High")

matriz_CD8
chi_matriz_CD8 <-chisq.test(matriz_CD8)

#CHI SQUARE ANALYSIS 
obsfreq <- c(2058,641)
nullprobs <- c(0.5,0.5)
#chisq.test(obsfreq,p=nullprobs)

obsfreq2 <- c(505,1141)
nullprobs <- c(0.5,0.5)

rows = 2 
matriz_CD8 = matrix(c(obsfreq,obsfreq2), 
                    nrow=rows,
                    byrow=TRUE)

rownames(matriz_CD8) = c("Untreated","ETP")
colnames(matriz_CD8) = c("Low", "High")

matriz_CD8
chi_matriz_CD8 <-chisq.test(matriz_CD8)

#GGPLOT 

labels <- m_CD8 %>% 
  group_by(treat) %>% 
  summarise(N=paste0("n =", n()))

geom_text(data = labels,aes(x=treat,y=5e+5,label=N), size=6) + 
  
  v_53BP1 <- ggplot(m_CD8, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  stat_summary(fun.data = n_fun, geom = "text") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("B lymphocytes") +
  annotate("text", x=1.45, y=50000, label = "23.7%", size = 6) +  #untreated and high
  annotate("text", x=1.45, y=10, label = "76.3%", size =6 ) + #untreated and low 
  annotate("text", x=2.43, y=50000, label = "69.3%", size =6 ) + #ETP and high 
  annotate("text", x=2.43, y=10, label = "30.7%", size = 6)



annotate("text", x=1.5, y=50000, label = "6.25%", size =6) +  #untreated and high
  annotate("text", x=1.5, y=2, label = "93.8%", size =6) + #untreated and low 
  annotate("text", x=2.43, y=50000, label = "51.9%", size =6) + #ETP and high 
  annotate("text", x=2.43, y=2, label = "48.1%", size =6)    



v_53BP11 <- ggplot(m_CD8, aes(x= treat,y=spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  geom_text(data = labels, aes(x=treat, y=5e+5, label=N), size= 6) +
  #stat_summary(fun.data = n_fun, geom = "text") +
  ylab("53BP1 Spot Intensity Sum log10(A.U) + 1") + 
  xlab("Treatment") + 
  #scale_fill_few(name ="Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("CD8+ T lymphocytes") +  
  annotate("text", x=1.45, y=50000, label = "23.7%", size = 6) +  #untreated and high
  annotate("text", x=1.45, y=10, label = "76.3%", size =6 ) + #untreated and low 
  annotate("text", x=2.43, y=50000, label = "69.3%", size =6 ) + #ETP and high 
  annotate("text", x=2.43, y=10, label = "30.7%", size = 6)  #EtP and low #


v_53BP11 <- ggplot(m_CD8, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  geom_text(data = labels,aes(x=treat, y=1e+6, label=N), size= 6) +
  #stat_summary(fun.data = n_fun, geom = "text") +
  #scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum log10(A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("CD8+ T lymphocytes") +
  labs(fill="Treatment") +
  annotate("text", x=1.42, y=50000, label = "35.7%") +  #untreated and high
  annotate("text", x=1.42, y=10, label = "64.3%") + #untreated and low 
  annotate("text", x=2.43, y=50000, label = "66.8%") + #ETP and high 
  annotate("text", x=2.43, y=10, label = "33.2%")  #EtP and low 


#STATS
#stats 
compare_means(spot_int53BP1_sum ~ treat,  data = High_spots_CD8, paired = FALSE)
my_comparisons = list( c("ETP", "Untreated") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(spot_int53BP1_sum ~ treat, data = High_spots_CD8)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(High_spots_CD8$spot_int53BP1_sum,    High_spots_CD8$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)


#print violin plot 
v_53BP11 + stat_compare_means(comparisons = my_comparisons, data = High_spots_CD8, paired = FALSE, label = "p.signif") + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + theme(axis.text = element_text(size=15))

#p.adj.signif
#+ geom_boxplot(width=0.1, notch = TRUE, boxlwd = 2)

#create violin plot for total number of spots for 53BP1
v_53BP1_spots <- ggplot(m_CD8, aes(x= treat, y= spots_53BP1, fill= treat)) +
  geom_violin() + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+2)) +
  stat_summary(fun.data = n_fun, geom = "text") +
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("CD8+ T lymphocytes") + 
  theme_bw()

#print violin plot 
v_53BP1_spots  + scale_fill_brewer(palette = "BuGn") + theme_bw()

# Box plot
#bp+scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest"))
# Scatter plot
#sp+scale_color_manual(values=wes_palette(n=3, name="GrandBudapest"))

#create biolin plot for yH2AX 
v_yH2AX <- ggplot(m_CD8, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  stat_summary(fun.data = n_fun, geom = "text") +
  ylab("yH2AX Nuclear Intensity Mean (A.U) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("CD8+ T lymphocytes")  

#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = m_CD8, paired = FALSE)
my_comparisons = list( c("ETP", "Untreated") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = m_CD8)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(m_CD8$nuc_intyH2AX_mean,    m_CD8$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")

#print violin plot with boxplot and significance values 
v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) 

#scatterplot 

scatter_CD8 <-ggplot(m_CD8, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) + 
  geom_point(shape=18, color="blue")+
  ggtitle("CD8") + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")




scatter_CD8 <-ggplot(m_CD8, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=treat)) + 
  geom_point(alpha = 0.5)+
  ggtitle("CD8") + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug()+
  scale_x_log10()+ 
  scale_y_log10()+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

scatter_B <-ggplot(m_B, aes(x= nuc_int53BP1_mean, y= nuc_intyH2AX_mean, color=treat)) + 
  geom_point(alpha = 0.5)+
  ggtitle("B-cells") + 
  xlab("53BP1 Mean Nuclear Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug()+
  scale_x_log10()+ 
  scale_y_log10()+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)


etp_CD8 <- filter(m_CD8, treat == "ETP")

etp_CD8_Reg <- lm(nuc_intyH2AX_mean ~ spot_int53BP1_sum, data = etp_CD8) 
summary(etp_CD8_Reg)

summary(etp_CD8_Reg)$r.squared



unt_CD8 <- filter(m_CD8, treat == "Untreated")

unt_CD8_Reg <- lm(nuc_intyH2AX_mean ~ spot_int53BP1_sum, data = unt_CD8) 
summary(unt_CD8_Reg)

summary(unt_CD8_Reg)$r.squared




#scatter_CD8 + facet_grid(vars(treat)) + labs(fill="Treatment")

scatter_CD8  + labs(fill="Treatment")
#----------------------------------------------
m_CD8_confilled <- ggplot(m_CD8, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat), show.legend = FALSE) +
  xlab("53BP1 sum spot intensity (a.u.)") + 
  ylab(expression(gamma ~ "H2AX nuclear intensity mean (a.u.)")) + 
  labs(fill="Treatment") + 
  ggtitle("CD8+T cells") + 
  scale_x_continuous(trans = 'log2',labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific)

#m_CD8_confilled +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
m_CD8_confilled +  theme_bw(base_size = 20) + theme(axis.text = element_text(size=15), legend.position = "bottom")    

fig3d_CD8 <- m_CD8_confilled + theme_classic(base_size = 20) + theme(axis.text = element_text(size = 15), legend.position = "bottom") 
