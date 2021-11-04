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
#set working directory 
#setwd("/Users/gallantkl")
#X190810_203 <- read_excel("Desktop/exp203_CD4_CD8_B_Monocytes_KG/data_203/190810-203-ALL-revised.xlsx")
#View(X190810_203)

X190810_203 <- read_excel("data/raw-data_v1/fig_3/190810-203-ALL-revised.xlsx")

#filter CD4 data 
m_CD4 <- filter(X190810_203, cell_type == "CD4") 
#, nuc_intyH2AX_mean >= 455.3

m_CD4$treat <- factor(m_CD4$treat, levels = c("Untreated","ETP"))

m_CD4_perc <- m_CD4 %>%
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

m_CD4_perc2 <- m_CD4 %>%
  group_by(treat, High_spots_yH2AX) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

perc_pos4 <-filter(m_CD4_perc, treat_lev == "High")

aggregate(spot_int53BP1_sum ~ treat, High_spots,mean)
aggregate(spot_int53BP1_sum ~ treat, High_spots,sd)
aggregate(nuc_intyH2AX_mean ~ treat, m_CD4,mean)
aggregate(nuc_intyH2AX_mean ~ treat, m_CD4,sd)
aggregate(nuc_intyH2AX_mean ~ treat,High_spots_yH2AX, mean)
table(High_spots_yH2AX$treat)
#print(m_CD4_perc)

#$res <- prop.test(x = c(933, 400), n = c(500, 500))
# Printing the results



bar_CD4 <- ggplot(m_CD4_perc, aes(x = treat, y = perc*100)) +
  geom_bar(
    aes(fill = factor(treat_lev, levels=c("Low","High"))),
    stat = "identity", position = position_stack()
  ) + 
  ggtitle("CD4") + 
  labs(x = "Treatment", y = "Total Cells (%)", fill = "53BP1 Spot Status") 
#geom_text(
#aes(label = perc*100, group = treat_lev)


bar_CD4 + scale_fill_brewer(palette = "BuGn") + theme_bw()


#create function to call sample sizes 
n_fun <- function(x){
  return(data.frame(y = max(x)*1.09, 
                    label = paste0("n = ",length(x))))
}

High_spots <- filter(m_CD4, treat_lev == "High")

High_spots_yH2AX <- filter(m_CD4,treat_yH2AX == "High")

High_CD4 <- High_spots %>%
  group_by(treat) %>% 
  get_summary_stats(spot_int53BP1_sum, type="median_mad")

yH2AX_CD4 <- m_CD4 %>%
  group_by(treat) %>% 
  get_summary_stats(nuc_intyH2AX_mean, type="median_mad")

#-------------------------------------------
#CHI SQUARE FOR yH2AX
#percentages 
j_CD4_perc <- m_CD4 %>%
  group_by(treat,treat_yH2AX ) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#CHI SQUARE ANALYSIS 
obsfreq <- c(2588,22)
nullprobs <- c(0.5,0.5)
#chisq.test(obsfreq,p=nullprobs)

obsfreq2 <- c(915,220)
nullprobs <- c(0.5,0.5)

rows = 2 
matriz_CD4 = matrix(c(obsfreq,obsfreq2), 
                    nrow=rows,
                    byrow=TRUE)

rownames(matriz_CD4) = c("DMSO","ETP")
colnames(matriz_CD4) = c("Low", "High")

matriz_CD4
chi_matriz_CD4 <-chisq.test(matriz_CD4)
#________________________________

#CHI SQUARE FOR 53BP1 
m_CD4perc <- m_CD4 %>%
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))
#CHI SQUARE ANALYSIS 
obsfreq <- c(1764,881)
nullprobs <- c(0.5,0.5)
#chisq.test(obsfreq,p=nullprobs)

obsfreq2 <- c(513,1141)
nullprobs <- c(0.5,0.5)

rows = 2 
matriz_CD4 = matrix(c(obsfreq,obsfreq2), 
                    nrow=rows,
                    byrow=TRUE)

rownames(matriz_CD4) = c("DMSO","ETP")
colnames(matriz_CD4) = c("Low", "High")

matriz_CD4
chi_matriz_CD4 <-chisq.test(matriz_CD4)


#CHI SQUARE ANALYSIS FOR CD4
#chi-square 
obsfreq <- c(1677,933)
nullprobs <- c(0.5,0.5)
chisq.test(obsfreq,p=nullprobs)

obsfreq2 <- c(377,758)
nullprobs <- c(0.5,0.5)

rows = 2 
matriz = matrix(c(obsfreq,obsfreq2), 
                nrow=rows,
                byrow=TRUE)

rownames(matriz) = c("Untreated","ETP")
colnames(matriz) = c("Low", "High")

matriz
chi_matriz <-chisq.test(matriz)
#chisq.test(obsfreq2,p=nullprobs)

labels_1 <-m_CD4 %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))

#geom_text(data = labels,aes(x=treat, y=5e+6, label=N), size= 6) +

#dodge <- position_dodge(width = 0.4)
#create violin plot for 53BP1 
v_53BP1_2 <- ggplot(m_CD4, aes(x= treat, y = spot_int53BP1_sum +1, fill = treat)) + 
  geom_violin(size=0.75, show.legend = FALSE) + 
  scale_y_continuous(trans= "log10", limits = c(1e+0, 1e+8)) + 
  #geom_text(data = labels_1, aes(x=treat, y=1e+8, label=N), size=6) + 
  geom_text(data = labels_1, aes(x=treat, y=5e+6, label=N), size=6) + 
  ylab("53BP1 spot intensity sum (a.u) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("CD4+T cells") + 
  annotate("text", x=1.40, y=50000, label = "35.7%", size =6) +  #untreated and high
  annotate("text", x=1.40, y=10, label = "64.3%", size =6) + #untreated and low 
  annotate("text", x=2.43, y=50000, label = "66.8%", size =6) + #ETP and high 
  annotate("text", x=2.43, y=10, label = "33.2%", size =6)  #EtP and low 


#theme(axis.text.x = element_text(size = 20)) 
#stats 
compare_means(spot_int53BP1_sum ~ treat,  data = High_spots, paired = FALSE)
my_comparisons = list( c("ETP", "Untreated") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(spot_int53BP1_sum ~ treat, data = High_spots)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(High_spots$spot_int53BP1_sum,    High_spots$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#print violin plot 
#v_53BP1_2  + stat_compare_means(comparisons = my_comparisons, label = "p.signif") + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + theme(axis.text = element_text(size=15)) 

fig3b_cd4 <- v_53BP1_2 + scale_fill_brewer(palette = "BuGn") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) +
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank())


#### SPOTS###
v_53BP1 <- ggplot(High_spots, aes(x= treat, y= spots_53BP1 +1, fill= treat)) +
  geom_violin() + 
  scale_y_continuous(trans = "log10", limits = c(1e0, 1e+1)) +
  stat_summary(fun.data = n_fun, geom = "text") +
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("CD4+ T lymphocytes") +
  labs(fill="Treatment") +
  theme(axis.line = element_line(colour="black", size = 3), panel.border = element_rect(colour = "black", fill=NA, size =2))

#print violin plot 
v_53BP1  + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20)


# Box plot
#bp+scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest"))
# Scatter plot
#sp+scale_color_manual(values=wes_palette(n=3, name="GrandBudapest"))

#create biolin plot for yH2AX 
v_yH2AX <- ggplot(m_CD4, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels_1, aes(x=treat, y=2e+5, label=N), size= 6) +
  #stat_summary(fun.data = n_fun, geom = "text") +
  labs(col = "Treatment") +
  ylab(expression(gamma ~ "H2AX nuclear intensity mean (a.u) + 1")) + 
  xlab("Treatment") + 
  labs(fill="Treatment") +
  ggtitle("CD4+T cells")  


#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = m_CD4, paired = FALSE)
my_comparisons = list( c("ETP", "Untreated") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = m_CD4)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(m_CD4$nuc_intyH2AX_mean,    m_CD4$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

untreat_CD4 <- subset(m_CD4, treat == "Untreated")
mean(untreat_CD4$nuc_intyH2AX_mean)
sd(untreat_CD4$nuc_intyH2AX_mean)
#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")

#print violin plot with boxplot and significance values 
#v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = 20) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + theme(axis.text= element_text(size=15)) 

fig3c_CD4 <- v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + 
  scale_fill_brewer(palette = "BuGn") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) +
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank())

#scatterplot formula
eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}

#scatterplot 
scatter_CD4 <-ggplot(m_CD4, aes(x= spot_int53BP1_sum, y= spot_intyH2AX_mean, color=treat)) + 
  geom_point(alpha = 0.5)+
  ggtitle("CD4") + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Spot Intensity") + 
  geom_rug() +
  scale_x_log10() + 
  scale_y_log10() + 
  geom_smooth(method= lm, se=FALSE, fullrange=TRUE)  

etp_CD4 <- filter(m_CD4, treat == "ETP")

etp_CD4_Reg <- lm(spot_intyH2AX_mean ~ spot_int53BP1_sum, data = etp_CD4) 
summary(etp_CD4_Reg)

summary(etp_CD4_Reg)$r.squared



unt_CD4 <- filter(m_CD4, treat == "Untreated")

unt_CD4_Reg <- lm(spot_intyH2AX_mean ~ spot_int53BP1_sum, data = unt_CD4) 
summary(unt_CD4_Reg)

summary(unt_CD4_Reg)$r.squared



#scatter_CD4 + facet_grid(vars(treat)) + labs(fill="Treatment") 

scatter_CD4 + labs(fill="Treatment") 

#----------------



#View(donor_untreated)

#created scatterplot with untreated data-> categorical 
m_CD4_all <- ggplot(m_CD4, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=cat_donor)) + 
  geom_point(alpha = 0.5)+
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug() +
  scale_x_continuous(trans = log2_trans()) + 
  scale_y_continuous(trans = log2_trans())

#run untreated data scatterplot 
m_CD4_all + facet_grid(treat ~ donor) + scale_color_grey(start = 0.3, end = 0.3) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  



#__________________________________________________________________________________________________

#geom_hexbin 

m_CD4_hex <- ggplot(m_CD4, mapping = aes(x= spot_int53BP1_sum , y= nuc_intyH2AX_mean)) + 
  geom_hex(bins = 60) +
  scale_fill_viridis_c() + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  scale_x_continuous(trans = 'log2') + 
  scale_y_continuous(trans = 'log2')

m_CD4_hex

#run geom hexbin
m_CD4_hex + facet_grid(treat ~ donor) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#--------
#geom_hexbar with both (facet_warap) 

#geom_hexbin 

m_CD4_hex <- ggplot(m_CD4, mapping = aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) + 
  geom_hex(bins = 20) +
  scale_fill_viridis_c() + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  scale_x_continuous(trans = log2_trans()) + 
  scale_y_continuous(trans = log2_trans())

m_CD4_hex

#run untreated data scatterplot 
m_CD4_hex + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  


#-----------
#contour plot 
#m_CD4$cat_donor<-cut(m_CD4$donor, seq(1,11,1), right=FALSE, labels=c(1:10))


#View(donor_untreated)

#created scatterplot with untreated data-> categorical 
m_CD4_contour <- ggplot(m_CD4, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color = treat)) + 
  geom_density2d(bins = 50) +
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  labs(fill="Treatment") + 
  #geom_rug() +
  scale_x_continuous(trans = 'log2') + 
  scale_y_continuous(trans = 'log2')

m_CD4_contour

#run untreated data scatterplot 
m_CD4_contour + facet_wrap(~donor, nrow = 2, ncol = 5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

#_______________________________
#another filled density/contour plot. 
m_CD4_confilled <- ggplot(m_CD4, aes(x= spot_int53BP1_sum , y= nuc_intyH2AX_mean )) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat), show.legend = FALSE) +
  xlab("53BP1 sum spot intensity (a.u.)") + 
  ylab(expression(gamma ~ "H2AX nuclear intensity mean (a.u)")) + 
  labs(fill="Treatment") + 
  ggtitle("CD4+T cells") + 
  scale_x_continuous(trans = 'log2', labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific)

#m_CD4_confilled +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
#m_CD4_confilled +  theme_bw(base_size = 20) + theme(axis.text = element_text(size=15))  
fig3d_cd4 <- m_CD4_confilled + theme_classic(base_size = 20) + theme(axis.text = element_text(size = 15), legend.position = "bottom") 

m_CD4_confilled <- ggplot(m_CD4, aes(x= spots_53BP1, y= nuc_intyH2AX_mean)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat)) +
  xlab("53BP1 Spots Per Cell") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  labs(fill="Treatment") + 
  ggtitle("CD4") + 
  scale_x_continuous(trans = 'log2') + 
  scale_y_continuous(trans = 'log2')

m_CD4_confilled +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  


#----- 
#yh2ax cd4

labels_yH2AX_CD4 <-High_spots_yH2AX %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))



v_yH2AX <- ggplot(High_spots_yH2AX, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels_yH2AX_CD4, aes(x=treat,y=5e+5,label=N), size=6) + 
  scale_fill_few(name="Treatment") +
  ylab("Pos-yH2AX Nuclear Intensity Mean log10(A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("CD4+ T lymphocytes") + 
  labs(fill="Treatment")  


#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = High_spots_yH2AX , paired = FALSE)
my_comparisons = list( c("Untreated", "ETP") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = High_spots_yH2AX)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(High_spots_yH2AX$nuc_intyH2AX_mean,    High_spots_yH2AX$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")

#print violin plot with boxplot and significance values 
v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text= element_text(size=15)) 


