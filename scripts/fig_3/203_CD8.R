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
X190810_203 <- read_excel("Desktop/exp203_CD4_CD8_B_Monocytes_KG/data_203/190810-203-ALL-revised.xlsx")
View(X190810_203)

#filter CD8 data 
m_CD8 <- filter(X190810_203, cell_type == "CD8") 


m_CD8$treat <- factor(m_CD8$treat, levels = c("Untreated","ETP"))

High_spots_yH2AX <- filter(m_CD8,treat_yH2AX == "High")

#percentage filter 
m_CD8_perc <- m_CD8 %>%
group_by(treat, treat_lev) %>%
  dplyr::summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

perc_pos5 <-filter(m_CD8_perc, treat_lev == "High")

aggregate(spot_int53BP1_sum ~ treat, m_CD8,mean)
aggregate(spot_int53BP1_sum ~ treat, m_CD8,sd)
aggregate(nuc_intyH2AX_mean ~ treat, m_CD8,mean)
aggregate(nuc_intyH2AX_mean ~ treat, m_CD8,sd)

table(High_spots_yH2AX$treat)

#barplot
bar_CD8 <- ggplot(m_CD8_perc, aes(x = treat, y = perc*100)) +
  geom_bar(
    aes(fill = factor(treat_lev, levels=c("Low","High"))),
    stat = "identity", position = position_stack()
  ) + 
  ggtitle("CD8+ T lymphocytes") + 
  labs(x = "Treatment", y = "Total Cells (%)", fill = "53BP1 Spot Status") 
#geom_text(
#aes(label = perc*100, group = treat_lev)



bar_CD8 + scale_fill_brewer(palette = "BuGn") + theme_bw()


#create function to call sample sizes 
n_fun <- function(x){
  return(data.frame(y = max(x)*1.1, 
                    label = paste0("n = ",length(x)), size = 6))
}

labels <-m_CD8 %>% 
  group_by(treat) %>% 
  summarise(N=paste0("n =", n()))




geom_text(data = labels,aes(x=treat,y=5e+5,label=N), size=6) + 

#create violin plot for 53BP1 
v_53BP1 <- ggplot(m_CD8, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(size=0.75) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  geom_text(data = labels,aes(x=treat,y=5e+5,label=N), size=6) + 
  #stat_summary(fun.data = n_fun, geom = "text") +
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum log10(A.U) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("CD8+ T lymphocytes") + 
  annotate("text", x=1.42, y=50000, label = "61.8%", size = 6) +  #untreated and high
  annotate("text", x=1.42, y=10, label = "38.2%", size = 6) + #untreated and low 
  annotate("text", x=2.43, y=50000, label = "71.2%", size = 6) + #ETP and high 
  annotate("text", x=2.43, y=10, label = "28.8%", size = 6)  #EtP and low 

  

#print violin plot 
v_53BP1 + scale_fill_brewer(palette = "BuGn") + theme_bw()

##create violin plot for high spots 53BP1
High_spots_CD8 <- filter(m_CD8, treat_lev == "High")

High_CD8 <- High_spots_CD8 %>%
  group_by(treat) %>% 
  get_summary_stats(spot_int53BP1_sum, type="median_mad")

yH2AX_CD8 <- m_CD8 %>%
  group_by(treat) %>% 
  get_summary_stats(nuc_intyH2AX_mean, type="median_mad")
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
  dplyr::summarise(N=paste0("n =", n()))

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
  
  

v_53BP1_cd8_ <- ggplot(m_CD8, aes(x= treat,y=spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  geom_text(data = labels, aes(x=treat, y=5e+6, label=N), size= 6) +
  #stat_summary(fun.data = n_fun, geom = "text") +
  ylab("53BP1 Spot Intensity Sum log10(A.U) + 1") + 
  xlab("Treatment") + 
  #scale_fill_few(name ="Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("CD8+ T lymphocytes") +  
  annotate("text", x=1.45, y=50000, label = "11.8%", size = 6) +  #untreated and high
  annotate("text", x=1.45, y=10, label = "88.2.%", size =6 ) + #untreated and low 
  annotate("text", x=2.43, y=50000, label = "52.6%", size =6 ) + #ETP and high 
  annotate("text", x=2.43, y=10, label = "47.4%", size = 6)  #EtP and low #


v_53BP11 <- ggplot(m_CD8, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e0, 1e+8)) +
  geom_text(data = labels,aes(x=treat, y=5e+6, label=N), size= 6) +
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
v_53BP1_cd8_ + stat_compare_means(comparisons = my_comparisons, label = "p.signif") + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + theme(axis.text = element_text(size=15))

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
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+5)) +
  stat_summary(fun.data = n_fun, geom = "text") +
  ylab("yH2AX Nuclear Intensity Mean log10(A.U) + 1") + 
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
  xlab("53BP1 Sum Spot Intensity log2(A.U.)") + 
  ylab("yH2AX Mean Nuclear Intensity log2(A.U.)") + 
  labs(fill="Treatment") + 
  ggtitle("CD8+ T lymphocytes") + 
  scale_x_continuous(trans = 'log2',labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific)

#m_CD8_confilled +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
m_CD8_confilled +  theme_bw(base_size = 20) + theme(axis.text = element_text(size=15), legend.position = "bottom")    

