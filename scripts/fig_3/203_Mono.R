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
View(X190810_203)


X190810_203 <- read_excel("data/raw-data_v1/fig_3/190810-203-ALL-revised.xlsx")

#filter Mono data 
m_Mono <- filter(X190810_203, cell_type == "Monocyte")

m_Mono$treat <- factor(m_Mono$treat, levels = c("Untreated","ETP"))

##SPOTS FOR YH@@AX) 
High_spots_yH2AX_m <- filter(m_Mono,treat_yH2AX == "High")
table(High_spots_yH2AX$treat)


m_Mono_sum <- m_Mono %>%
  group_by(treat) %>% 
  get_summary_stats(spot_int53BP1_sum, type="mean_sd")


#percentage 
m_Mono_perc <- m_Mono %>%
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

perc_posM <-filter(m_Mono_perc, treat_lev == "High")

aggregate(spot_int53BP1_sum ~ treat, High_spots_mono,mean)
aggregate(spot_int53BP1_sum ~ treat, High_spots_mono,sd)
aggregate(nuc_intyH2AX_mean ~ treat, m_Mono,mean)
aggregate(nuc_intyH2AX_mean ~ treat, m_Mono,sd)


table(High_spots_yH2AX$treat)
b#create function to call sample sizes 
n_fun <- function(x){
  return(data.frame(y = max(x)*1.1, 
                    label = paste0("n = ",length(x))))
}

#create violin plot for 53BP1 
v_53BP1 <- ggplot(m_Mono, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin() + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  stat_summary(fun.data = n_fun, geom = "text") +
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("Monocyte") + 
  theme_bw()

#print violin plot 
v_53BP1 + scale_fill_brewer(palette = "BuGn") + theme_bw()

High_spots_mono <- filter(m_Mono, treat_lev == "High")

#calculate medians 
High_Mono <- High_spots_mono %>%
  group_by(treat) %>% 
  get_summary_stats(spot_int53BP1_sum, type="median_mad")


yH2AX_Mono <- m_Mono %>%
  group_by(treat) %>% 
  get_summary_stats(nuc_intyH2AX_mean, type="median_mad")

#CHI SQUARE FOR yH2AX
#percentages 
j_Mono_perc <- m_Mono %>%
  group_by(treat,treat_yH2AX ) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#CHI SQUARE ANALYSIS 
obsfreq <- c(494,55)
nullprobs <- c(0.5,0.5)
#chisq.test(obsfreq,p=nullprobs)

obsfreq2 <- c(42,373)
nullprobs <- c(0.5,0.5)

rows = 2 
matriz_Mono = matrix(c(obsfreq,obsfreq2), 
                    nrow=rows,
                    byrow=TRUE)

rownames(matriz_Mono) = c("DMSO","ETP")
colnames(matriz_Mono) = c("Low", "High")

matriz_Mono
chi_matriz_Mono <-chisq.test(matriz_Mono)

#CHI SQUARE ANALYSIS  (low/high)
obsfreq <- c(516,33)
nullprobs <- c(0.5,0.5)
#chisq.test(obsfreq,p=nullprobs)

obsfreq2 <- c(376,39)
nullprobs <- c(0.5,0.5)

rows = 2 
matriz_Mono = matrix(c(obsfreq,obsfreq2), 
                  nrow=rows,
                  byrow=TRUE)

rownames(matriz_Mono) = c("Untreated","ETP")
colnames(matriz_Mono) = c("Low", "High")

matriz_Mono
chi_matriz_Mono <-chisq.test(matriz_Mono)
chi_matriz_Mono
chi_matriz_Mono$p.value

labels <-m_Mono %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))



#create violin plot for mean nuclear intensity for 53BP1 
v_53BP1_2 <- ggplot(m_Mono, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  geom_text(data = labels,aes(x=treat, y=5e+6, label=N), size= 6) +
  #stat_summary(fun.data = n_fun, geom = "text") +
  ylab("53BP1 spot intensity sum (a.u.) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("Monocytes") + 
  annotate("text", x=1.25, y=50000, label = "6.0%", size =6) +  #untreated and high
  annotate("text", x=1.25, y=10, label = "94.0%", size =6 ) + #untreated and low 
  annotate("text", x=2.25, y=50000, label = "9.4%", size =6) + #ETP and high 
  annotate("text", x=2.25, y=10, label = "90.6%", size =6)  #EtP and low 
 
compare_means(spot_int53BP1_sum ~ treat,  data = High_spots_mono, paired = FALSE)
my_comparisons = list( c("ETP", "Untreated") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(spot_int53BP1_sum ~ treat, data = High_spots_mono)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(High_spots_mono$spot_int53BP1_sum,    High_spots_mono$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)
#print violin plot with stats 
v_53BP1_2  + stat_compare_means(comparisons = my_comparisons, data= High_spots_mono,  label = "p.signif") + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + theme(axis.text = element_text(size=15))
#+ geom_boxplot(width=0.1, notch = TRUE, boxlwd = 2)

fig3b_mono <- v_53BP1_2 + scale_fill_brewer(palette = "BuGn") + 
  stat_compare_means(comparisons = my_comparisons, data = High_spots_mono, label = "p.signif") + 
  theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) +
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank())
#print violin plot 
#v_53BP1_2 + geom_boxplot(width=0.1) + scale_fill_brewer(palette = "Set2")

#boxplot 
bar_Mono <- ggplot(m_Mono_perc, aes(x = treat, y = perc*100)) +
  geom_bar(
    aes(fill = factor(treat_lev, levels=c("Low","High"))),
    stat = "identity", position = position_stack()
  ) + 
  ggtitle("Monocyte") + 
  labs(x = "Treatment", y = "Percentage of Total Cells (%)", fill = "53BP1 Spot Status") 


bar_Mono + scale_fill_brewer(palette = "BuGn") + theme_bw()
  #geom_text(
    #aes(label = perc*100, group = treat_lev)
  )

# Box plot
#bp+scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest"))
# Scatter plot
#sp+scale_color_manual(values=wes_palette(n=3, name="GrandBudapest"))

#create violin plot for yH2AX 
v_yH2AX <- ggplot(m_Mono, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels,aes(x=treat, y=2e+5, label=N), size= 6) +
  #stat_summary(fun.data = n_fun, geom = "text") +
  labs(col = "Treatment") +
  ylab(expression(gamma ~ "H2AX nuclear intensity mean (a.u.) + 1")) + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("Monocytes")  


#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = m_Mono, paired = FALSE)
my_comparisons = list( c("ETP", "Untreated") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = m_Mono)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(m_Mono$nuc_intyH2AX_mean,    m_Mono$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")

#print violin plot with boxplot and significance values 
v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + theme(axis.text=element_text(size=15)) 

fig3c_mono <- v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + scale_fill_brewer(palette = "BuGn") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) +
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank())


#scatterplot 
scatter_mono <-ggplot(m_Mono, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=treat)) + 
  geom_point(alpha = 0.5)+
  ggtitle("Monocytes") + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug()+
  scale_x_log10()+ 
  scale_y_log10()+ 
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")



#scatter_mono + facet_grid(vars(treat)) + labs(fill="Treatment") 
scatter_mono + labs(fill="Treatment") 

#--------- 

#yH2AX violin plot redo 

labels_yH2AX_m <-High_spots_yH2AX_m %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))



v_yH2AX <- ggplot(High_spots_yH2AX_m, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels_yH2AX_m, aes(x=treat,y=2e+5,label=N), size=6) + 
  scale_fill_few(name="Treatment") +
  ylab("Pos-yH2AX Nuclear Intensity Mean log10(A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("Monocytes") + 
  labs(fill="Treatment")  


#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = High_spots_yH2AX_m, paired = FALSE)
my_comparisons = list( c("Untreated", "ETP") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = High_spots_yH2AX_m)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(High_spots_yH2AX_m$nuc_intyH2AX_mean,    High_spots_yH2AX_m$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")

#print violin plot with boxplot and significance values 
v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text= element_text(size=15)) 

#-----------

m_mono_confilled <- ggplot(m_Mono, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat), show.legend = FALSE) +
  xlab("53BP1 sum spot intensity (a.u.)") + 
  ylab(expression(gamma ~ "H2AX nuclear intensity mean (a.u.)")) + 
  labs(fill="Treatment") + 
  ggtitle("Monocytes") + 
  scale_x_continuous(trans = 'log2',labels =scales::scientific, limits = c(2e+4, 2e+5)) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific)

#m_CD8_confilled +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
#m_CD8_confilled +  theme_bw(base_size = 20) + theme(axis.text = element_text(size=15), legend.position = "bottom")    

fig3d_mono <- m_mono_confilled + theme_classic(base_size = 20) + theme(axis.text = element_text(size = 15), legend.position = "bottom") 

