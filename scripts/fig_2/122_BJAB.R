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
setwd("/Users/gallantkl")

X_119_BJAB <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190124-119-KLG-Jurkat-BJAB-THP1-Ad2x_20190124_113931[3510]/190218-119-BJAB.xlsx")
View(X_119_BJAB)

#library(readxl)
#x_119_BJAB <- read_excel("Desktop/190218-121-122-BJAB-Jurkat-ETPtitration-BJAB/190218-122-KG-BJAB-yH2AX-53bp1-dose_20190218_150512[3608]/190218-119-BJAB.xlsx")
#View(x_119_BJAB)
#filter 0 and 30 etp  
j_fil_BJAB <- X_119_BJAB %>%
  filter(treat %in% c("DMSO", "30"))

#reorder groups
j_fil_BJAB$treat <- factor(j_fil_BJAB$treat, levels = c("DMSO","30"))

#find counts and percentages 
j_fil_BJABperc <- j_fil_BJAB %>%
  group_by(treat, treat_lev) %>%
  dplyr::summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

perc_pos2 <- filter(j_fil_BJABperc, treat_lev == "Positive")



high_spots_yH2AX2 <- filter(j_fil_BJAB, treat_yH2AX == "High")

table(high_spots_yH2AX2$treat)

#barplot
b_barplot <- ggplot(j_fil_BJABperc, aes(x = treat, y = perc*100)) +
  geom_bar(
    aes(fill = treat_lev),
    stat = "identity", position = position_stack()
  ) + 
  ggtitle("BJAB") + 
  labs(x = "Treatment", y = "Total Cells (%)", fill = "53BP1 Spot Status")  

b_barplot + scale_fill_brewer(palette = "BuGn") + theme_bw()

rperc_pos <- filter(j_fil_BJABperc, treat_lev == "Positive")

j_fil_BJABperc

#------------------------------------------

#Chi-square for yH2AX

#percentages 
j_BJAB_perc <- j_fil_BJAB %>%
  group_by(treat,treat_yH2AX ) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#CHI SQUARE ANALYSIS 
obsfreq <- c(5183,314)
nullprobs <- c(0.5,0.5)
#chisq.test(obsfreq,p=nullprobs)

obsfreq2 <- c(1492,4701)
nullprobs <- c(0.5,0.5)

rows = 2 
matriz_BJAB = matrix(c(obsfreq,obsfreq2), 
                     nrow=rows,
                     byrow=TRUE)

rownames(matriz_BJAB) = c("DMSO","ETP")
colnames(matriz_BJAB) = c("Low", "High")

matriz_BJAB
chi_matriz_BJAB <-chisq.test(matriz_BJAB)
#------------

#CHI SQUARE FOR 53BP1 
j_fil_BJABperc <- j_fil_BJAB %>%
  group_by(treat, treat_lev) %>%
  dplyr::summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))
#CHI SQUARE ANALYSIS 
obsfreq <- c(5197,300)
nullprobs <- c(0.5,0.5)
#chisq.test(obsfreq,p=nullprobs)

obsfreq2 <- c(3603,2590)
nullprobs <- c(0.5,0.5)

rows = 2 
matriz_BJAB = matrix(c(obsfreq,obsfreq2), 
                    nrow=rows,
                    byrow=TRUE)

rownames(matriz_BJAB) = c("DMSO","ETP")
colnames(matriz_BJAB) = c("Low", "High")

matriz_BJAB
chi_matriz_BJAB <-chisq.test(matriz_BJAB)


#--------------------
#filled contour plots 

BJAB_confilled <- ggplot(j_fil_BJAB, aes(x= spot_int53BP1_sum , y= nuc_intyH2AX_mean)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat), show.legend = FALSE) +
  xlab("53BP1 Sum Spot Intensity log2(A.U.)") + 
  ylab("yH2AX Mean Nuclear Intensity log2(A.U.)") + 
  labs(fill="Treatment") + 
  ggtitle("BJAB") + 
  scale_x_continuous(trans = 'log2', labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific)

#BJAB_confilled +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
BJAB_confilled + theme_bw(base_size = 20) + theme(axis.text = element_text(size=15)) 

View(j_fil_BJAB)
#names <- c("DMSO", "30")

#lev1 <- c("DMSO", "30")

#j_fil_cat <- j_fil %>%
#mutate(names = factor(names, levels = lev1))


#msleep %>% 
#select(order, name, sleep_total) %>% 
#filter(order %in% c("Didelphimorphia", "Diprotodontia"))
#df <- tibble(names, values) %>%
#mutate(names = factor(names, levels = lev1))

#create function to call sample sizes 
n_fun <- function(x){
  return(data.frame(y = max(x)*1.1, 
                    label = paste0("n = ",length(x)), size = 6))
}

#stats
High_spots_BJAB <- filter(j_fil_BJAB, treat_lev == "Positive")

High_BJAB <- High_spots_BJAB %>%
  group_by(treat) %>% 
  get_summary_stats(spot_int53BP1_sum, type="median_mad")

yH2AX_BJAB <- j_fil_BJAB %>%
  group_by(treat) %>% 
  get_summary_stats(nuc_intyH2AX_mean, type="median_mad")

labels3 <-j_fil_BJAB %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))


#create violin plot for 53BP1 
v_53BP1 <- ggplot(j_fil_BJAB, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+9)) +
  geom_text(data = labels,aes(x=treat,y=1e+8,label=N), size=6) + 
  #stat_summary(fun.data = n_fun, geom = "text") +
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum log10(A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("BJAB") +
  labs(fill="Treatment") + 
  annotate("text", x=1.35, y=50000, label = "5.5%", size = 6) +  #untreated and high
  annotate("text", x=1.35, y=10, label = "94.5%", size = 6) + #untreated and low 
  annotate("text", x=2.35, y=50000, label = "41.8%", size = 6) + #ETP and high 
  annotate("text", x=2.35, y=10, label = "58.2%", size = 6)  #EtP and low 

#stats 
compare_means(spot_int53BP1_sum ~ treat,  data = High_spots_BJAB, paired = FALSE)
my_comparisons = list( c("30", "DMSO") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(spot_int53BP1_sum ~ treat, data = High_spots_BJAB)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(High_spots_BJAB$spot_int53BP1_sum,    High_spots_BJAB$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)
#print violin plot with stats 
v_53BP1  + stat_compare_means(comparisons = my_comparisons, label = "p.signif") + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text = element_text(size=15))


#+ geom_boxplot(width=0.1, notch = TRUE, boxlwd = 2)


#create violin plot for 53BP1 
v_53BP1_B <- ggplot(j_fil_BJAB, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin() + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  stat_summary(fun.data = n_fun, geom = "text") +
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  xlab("Treatment") + 
  theme_bw()

#stats 

#print violin plot 
v_53BP1_B + scale_fill_brewer(palette = "Set2")


aggregate(spot_int53BP1_sum ~ treat, j_fil_BJAB,mean)
aggregate(spot_int53BP1_sum ~ treat, j_fil_BJAB,sd)

aggregate(spot_int53BP1_sum ~ treat, High_spots_BJAB,mean)
aggregate(spot_int53BP1_sum ~ treat, High_spots_BJAB,sd)
aggregate(nuc_intyH2AX_mean ~ treat, j_fil_BJAB,mean)
aggregate(nuc_intyH2AX_mean ~ treat, j_fil_BJAB,sd)
# Box plot
#bp+scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest"))
# Scatter plot
#sp+scale_color_manual(values=wes_palette(n=3, name="GrandBudapest"))

#create biolin plot for yH2AX 
v_yH2AX <- ggplot(j_fil_BJAB, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels3,aes(x=treat,y=1e+5,label=N), size=6) + 
  #stat_summary(fun.data = n_fun, geom = "text") +
  scale_fill_few(name="Treatment") +
  ylab("yH2AX Nuclear Intensity Mean log10(A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("BJAB") + 
  labs(fill="Treatment")  


#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = j_fil_BJAB, paired = FALSE)
my_comparisons = list( c("DMSO", "30") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = j_fil_BJAB)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(j_fil_BJAB$nuc_intyH2AX_mean,    j_fil_BJAB$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")

#print violin plot with boxplot and significance values 
v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text= element_text(size=15)) 

scatter_BJAB <-ggplot(j_fil_BJAB, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) + 
  geom_point(shape=18, color="blue")+
  ggtitle("BJAB") + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")


scatter_BJAB<-ggplot(j_fil_BJAB, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=treat)) + 
  geom_point(alpha = 0.5)+
  ggtitle("BJAB") + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug() +
  scale_x_log10() + 
  scale_y_log10() + 
  geom_smooth(method= lm, se=FALSE, fullrange=TRUE)  

etp_BJAB <- filter(j_fil_BJAB, treat == "30")

etp_BJAB_Reg <- lm(nuc_intyH2AX_mean ~ spot_int53BP1_sum, data = etp_BJAB) 
summary(etp_BJAB_Reg)

summary(etp_BJAB_Reg)$r.squared



unt_BJAB <- filter(j_fil_BJAB, treat == "DMSO")

unt_BJAB_Reg <- lm(nuc_intyH2AX_mean ~ spot_int53BP1_sum, data = unt_BJAB) 
summary(unt_BJAB_Reg)

summary(unt_BJAB_Reg)$r.squared



#scatter_BJAB + facet_grid(vars(treat)) + labs(fill="Treatment")

scatter_BJAB + labs(fill="Treatment")


#------
#yH2AX redo 

labels_yH2AX_BJAB <-high_spots_yH2AX2 %>% 
  group_by(treat) %>% 
  summarise(N=paste0("n =", n()))

geom_text(data = labels_yH2AX_BJAB,aes(x=treat,y=5e+4,label=N), size=6) + 

v_yH2AX1 <- ggplot(high_spots_yH2AX2, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+5)) +
  geom_text(data = labels_yH2AX_BJAB,aes(x=treat,y=5e+4,label=N), size=6) + 
  scale_fill_few(name="Treatment") +
  ylab("Pos-yH2AX Nuclear Intensity Mean log10(A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("BJAB") + 
  labs(fill="Treatment")  


#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = high_spots_yH2AX2, paired = FALSE)
my_comparisons = list( c("DMSO", "30") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = high_spots_yH2AX2)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(high_spots_yH2AX2$nuc_intyH2AX_mean,    high_spots_yH2AX2$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")

#print violin plot with boxplot and significance values 
v_yH2AX1 + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text= element_text(size=15)) 
