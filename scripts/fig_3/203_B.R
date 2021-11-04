#the purpose of this code is for data visualization of cells 
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
#View(X190810_203)

X190810_203 <- read_excel("data/raw-data_v1/fig_3/190810-203-ALL-revised.xlsx")

#filter B data 
m_B <- filter(X190810_203, cell_type == "B-cell")

m_B$treat <- factor(m_B$treat, levels = c("Untreated","ETP"))


High_spots_yH2AX_ <- filter(m_B,treat_yH2AX == "High")

#percentage 
m_B_perc <- m_B %>% 
group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

perc_pos7 <-filter(m_B_perc, treat_lev == "High")

aggregate(spot_int53BP1_sum ~ treat, High_spots_B,mean)
aggregate(spot_int53BP1_sum ~ treat, High_spots_B,sd)
aggregate(nuc_intyH2AX_mean ~ treat, m_B,mean)
aggregate(nuc_intyH2AX_mean ~ treat, m_B,sd)

table(High_spots_yH2AX$treat)


#create function to call sample sizes 
n_fun <- function(x){
  return(data.frame(y = max(x)*1.1, 
                    label = paste0("n = ",length(x))))
}

#CHI SQUARE ANALYSIS 


#create violin plot for 53BP1 
v_53BP1 <- ggplot(m_B, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  stat_summary(fun.data = n_fun, geom = "text") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("B cells") + 
  labs(fill="Treatment") + 
  annotate("text", x=1, y=50000, label = "6.25%", size =6) +  #untreated and high
  annotate("text", x=1, y=2, label = "93.8%", size =6 ) + #untreated and low 
  annotate("text", x=2, y=50000, label = "51.9%", size =6) + #ETP and high 
  annotate("text", x=2, y=2, label = "48.1%", size =6) #EtP and low 
  

#print violin plot 
v_53BP1 + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) 


#barplot B cell 
bar_B <- ggplot(m_B_perc, aes(x = treat, y = perc*100)) +
  geom_bar(
    aes(fill = factor(treat_lev, levels=c("Low","High"))),
    stat = "identity", position = position_stack()
  ) + 
  ggtitle("B lymphocytes") + 
  labs(x = "Treatment", y = "Percentage of Total Cells (%)", fill = "53BP1 Spot Status") 
#geom_text(
#aes(label = perc*100, group = treat_lev)

bar_B + scale_fill_brewer(palette = "BuGn") + theme_bw()



# Box plot
#bp+scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest"))
# Scatter plot
#sp+scale_color_manual(values=wes_palette(n=3, name="GrandBudapest"))
r <- m_B %>%
  group_by(treat) %>% 
  get_summary_stats(spot_int53BP1_sum, type="median_mad")

#filter out high intensity spots from B cells 
High_spots_B <- filter(m_B, treat_lev == "High")

#stats for each group

High_B <- High_spots_B %>%
  group_by(treat) %>% 
  get_summary_stats(spot_int53BP1_sum, type="median_mad")

yH2AX_B <- m_B %>%
  group_by(treat) %>% 
  get_summary_stats(nuc_intyH2AX_mean, type="median_mad")

#CHI SQUARE FOR yH2AX
#percentages 
j_B_perc <- m_B %>%
  group_by(treat,treat_yH2AX ) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#CHI SQUARE ANALYSIS 
obsfreq <- c(2480,12)
nullprobs <- c(0.5,0.5)
#chisq.test(obsfreq,p=nullprobs)

obsfreq2 <- c(1065,323)
nullprobs <- c(0.5,0.5)

rows = 2 
matriz_B = matrix(c(obsfreq,obsfreq2), 
                    nrow=rows,
                    byrow=TRUE)

rownames(matriz_B) = c("DMSO","ETP")
colnames(matriz_B) = c("Low", "High")

matriz_B
chi_matriz_B <-chisq.test(matriz_B)

#CHI SQUARE ANALYSIS  (low/high)
obsfreq <- c(2280,152)
nullprobs <- c(0.5,0.5)
#chisq.test(obsfreq,p=nullprobs)

obsfreq2 <- c(668,720)
nullprobs <- c(0.5,0.5)

rows = 2 
matriz_B = matrix(c(obsfreq,obsfreq2), 
                  nrow=rows,
                  byrow=TRUE)

rownames(matriz_B) = c("Untreated","ETP")
colnames(matriz_B) = c("Low", "High")

matriz_B
chi_matriz_B <-chisq.test(matriz_B)
chi_matriz_B
chi_matriz_B$p.value

labels <-m_B %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))

  
# use this one! 
  v_53BP1 <- ggplot(m_B, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  geom_text(data = labels,aes(x=treat, y=5e+6, label=N), size= 6) +
  #stat_summary(fun.data = n_fun, geom = "text") +
  #scale_fill_few(name="Treatment") +
  ylab("53BP1 spot intensity sum (a.u.) + 1") + 
  xlab("Treatment") + 
  ggtitle("B cells") +
  labs(fill="Treatment") +
  annotate("text", x=1.32, y=50000, label = "6.25%", size =6) +  #untreated and high
  annotate("text", x=1.32, y=10, label = "93.8%", size =6) + #untreated and low 
  annotate("text", x=2.32, y=50000, label = "51.9%", size =6) + #ETP and high 
  annotate("text", x=2.32, y=10, label = "48.1%", size =6)  #EtP and low   
  
  
  
  
  

v_53BP1 <- ggplot(m_B, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  stat_summary(fun.data = n_fun, geom = "text") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("B lymphocytes") +
  annotate("text", x=1.5, y=50000, label = "6.25%", size =6) +  #untreated and high
  annotate("text", x=1.5, y=2, label = "93.8%", size =6) + #untreated and low 
  annotate("text", x=2.43, y=50000, label = "51.9%", size =6) + #ETP and high 
  annotate("text", x=2.43, y=2, label = "48.1%", size =6)  
 #EtP and low 
 

#STATS
#stats 
compare_means(spot_int53BP1_sum ~ treat,  data = High_spots_B, paired = FALSE)
my_comparisons = list( c("ETP", "Untreated") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(spot_int53BP1_sum ~ treat, data = High_spots_B)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(High_spots_B$spot_int53BP1_sum,    High_spots_B$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)


#print violin plot 
v_53BP1 + stat_compare_means(comparisons = my_comparisons, label = "p.signif") + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + theme(axis.text = element_text(size=15))
#+ geom_boxplot(width=0.1, notch = TRUE, boxlwd = 2)
#medians
fig3b_bcell <- v_53BP1 + scale_fill_brewer(palette = "BuGn") + 
  stat_compare_means(comparisons = my_comparisons, data = High_spots_B, label = "p.signif") + 
  theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) +
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank())


#create violin plot for yH2AX 
v_yH2AX <- ggplot(m_B, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels,aes(x=treat, y=2e+5, label=N), size= 6) +
  #stat_summary(fun.data = n_fun, geom = "text") +
  labs(col = "Treatment") +
  ylab(gamma ~ "H2AX nuclear intensity mean (a.u) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("B cells") 

#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = m_B, paired = FALSE)
my_comparisons = list( c("ETP", "Untreated") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = m_B)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(m_B$nuc_intyH2AX_mean,    m_B$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")

#print violin plot with boxplot and significance values 
v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + theme(axis.text = element_text(size=15)) 

fig3c_bcell <- v_yH2AX + scale_fill_brewer(palette = "BuGn") + 
  geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) +
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank())

#scatterplot

scatter_B <-ggplot(m_B, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) + 
  geom_point(shape=18, color="blue")+
  ggtitle("B lymphocytes") + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")


scatter_B <-ggplot(m_B, aes(x= nuc_int53BP1_mean, y= nuc_intyH2AX_mean, color=treat)) + 
  geom_point(alpha = 0.5)+
  ggtitle("B lymphocytes") + 
  xlab("53BP1 Mean Nuclear Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug()+
  scale_x_log10()+ 
  scale_y_log10()+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)

etp_B <- filter(m_B, treat == "ETP")

etp_B_Reg <- lm(nuc_intyH2AX_mean ~ spot_int53BP1_sum, data = etp_B) 
summary(etp_B_Reg)

summary(etp_B_Reg)$r.squared



unt_B <- filter(m_B, treat == "Untreated")

unt_B_Reg <- lm(nuc_intyH2AX_mean ~ spot_int53BP1_sum, data = unt_B) 
summary(unt_B_Reg)

summary(unt_B_Reg)$r.squared

#scatter_B + facet_grid(vars(treat)) + labs(fill="Treatment") 

scatter_B  + labs(fill="Treatment") 

#----------------------------------------------
m_B_confilled <- ggplot(m_B, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat), show.legend = FALSE) +
  xlab("53BP1 sum spot intensity (a.u.)") + 
  ylab(gamma ~ "H2AX nuclear intensity mean (a.u)") + 
  labs(fill="Treatment") + 
  ggtitle("B cells") + 
  scale_x_continuous(trans = 'log2', labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2',labels =scales::scientific)

#m_B_confilled +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
#m_B_confilled + theme_bw(base_size = 20) + theme(axis.text = element_text(size=15))

fig3d_bcell <- m_B_confilled + theme_classic(base_size = 20) + theme(axis.text = element_text(size = 15), legend.position = "bottom") 
#----- 

#yH2AX violin plot redo 

labels_yH2AX_B <-High_spots_yH2AX_ %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))



v_yH2AX <- ggplot(High_spots_yH2AX_, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels_yH2AX_B, aes(x=treat,y=2e+5,label=N), size=6) + 
  scale_fill_few(name="Treatment") +
  ylab("Pos-yH2AX Nuclear Intensity Mean log10(A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("B lymphocytes") + 
  labs(fill="Treatment")  


#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = High_spots_yH2AX_, paired = FALSE)
my_comparisons = list( c("Untreated", "ETP") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = High_spots_yH2AX_)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(High_spots_yH2AX_$nuc_intyH2AX_mean,    High_spots_yH2AX_$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")

#print violin plot with boxplot and significance values 
v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text= element_text(size=15)) 




