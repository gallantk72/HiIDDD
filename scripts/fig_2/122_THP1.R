#the purpose of this code is for data visualization of cells 
library(tidyverse)
library(ggthemes)
library(scales)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(stats)
library(forcats)

#find and read data from the directory
library(readxl)
setwd("/Users/gallantkl")

X_122_THP1 <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190218-122-KG-thp1-yH2AX-53bp1-dose_20190218_150512[3608]/190218-EXP122-DoseResponse-THP1.xlsx")
View(X_122_THP1)

#filter 0 and 30 etp  
j_fil_THP <- X_122_THP1 %>%
  filter(treat %in% c("DMSO", "30"))

j_fil_THP$treat


recode(j_fil_THP$treat, DMSO = "DMSO", 30 = "ETP" )

levels(j_fil_THP$treat) <- c("DMSO", "ETP")

j_fil_THP$treat <- recode_factor(j_fil_THP$treat, DMSO = "DMSO", ETP = "30")

#reorder
j_fil_THP$treat <- factor(j_fil_THP$treat, levels = c("DMSO","30"))

#percentage
j_fil_THPperc <- j_fil_THP %>%
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

j_fil_THPperc2 <- j_fil_THP %>%
  group_by(treat, treat_yH2AX) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

-------

perc_pos3 <- filter(j_fil_THPperc, treat_lev == "High")
#barplot B cell 
t_bar <- ggplot(j_fil_THPperc, aes(x = treat, y = perc*100)) +
  geom_bar(aes(fill = treat_lev),
    stat = "identity", position = position_stack(reverse = TRUE), show.legend = TRUE) + 
  ggtitle("THP-1") +
  labs(x = "Treatment", y = "Percentage of Total Cells (%)", fill = "53BP1 Spot Status")  
#stat_bin(aes(label = sprintf("%.02f %%", perc*100), group = treat_lev), geom ="text")

#t_bar <- arrange(j_fil_THPperc, treat_lev, desc(perc*100))

t_bar + scale_fill_manual(values = c("#99d8c8","#e5f5f8")) + theme_bw(base_size = 20) + theme(axis.text = element_text(size=15)) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP"))  

#names <- c("DMSO", "30")

#lev1 <- c("DMSO", "30")

#j_fil_cat <- j_fil %>%
  mutate(names = factor(names, levels = lev1))


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

#high spots 
High_spots_THP <- filter(j_fil_THP, treat_lev == "High")

High_spots_yH2AX_THP <- filter(j_fil_THP, treat_yH2AX == "High")

#stat summary

aggregate(spot_int53BP1_sum +1 ~ treat, j_fil_THP,mean)
aggregate(spot_int53BP1_sum +1 ~ treat, j_fil_THP,sd)
aggregate(nuc_intyH2AX_mean +1  ~ treat, j_fil_THP,mean)
aggregate(nuc_intyH2AX_mean +1 ~ treat, j_fil_THP,sd)


High_THP <- High_spots_THP %>%
  group_by(treat) %>% 
  get_summary_stats(spot_int53BP1_sum, type="median_mad")

yH2AX_THP <- j_fil_THP %>%
  group_by(treat) %>% 
  get_summary_stats(nuc_intyH2AX_mean, type="median_mad")

labels <-j_fil_THP %>% 
  group_by(treat) %>% 
  summarise(N=paste0("n =", n()))
geom_text(data = labels,aes(x=treat,y=5e+5,label=N), size=6) + 

#create violin plot for 53BP1 
v_53BP1 <- ggplot(j_fil_THP, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  geom_text(data = labels,aes(x=treat,y=5e+6,label=N), size=6) + 
  #stat_summary(fun.data = n_fun, geom = "text") +
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum log10(A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("THP-1") +
  labs(fill="Treatment") + 
  annotate("text", x=1.39, y=80000, label = "61.8%", size = 6) +  #untreated and high
  annotate("text", x=1.39, y=10, label = "38.2%", size = 6) + #untreated and low 
  annotate("text", x=2.42, y=80000, label = "71.2%", size = 6) + #ETP and high 
  annotate("text", x=2.42, y=10, label = "28.8%", size = 6)  #EtP and low 
#stats 
compare_means(spot_int53BP1_sum ~ treat,  data = High_spots_THP, paired = FALSE)
my_comparisons = list( c("30", "DMSO") )
stat_compare_means(comparisons = my_comparisons)
c(mean = mean(x), sd = sd(x))
#results of means
aggregate(spot_int53BP1_sum ~ treat, High_spots_THP,mean)
aggregate(spot_int53BP1_sum ~ treat, High_spots_THP,sd)
#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(spot_int53BP1_sum ~ treat, data = High_spots_THP)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(High_spots_THP$spot_int53BP1_sum,   High_spots_THP$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)
#scale_x_discrete(labels=c("0.5" = "Dose 0.5", "1" = "Dose 1",
                          #"2" = "Dose 2"))
# look at the comparisons and p-values
head(tidy_ad_pairwise)
#print violin plot with stats 
v_53BP1 + stat_compare_means(comparisons = my_comparisons, label = "p.signif") + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text = element_text(size=15))
#print violin plot  
v_53BP1 + scale_fill_brewer(palette = "Set2")
geom_boxplot(width=0.1, notch = TRUE, boxlwd = 2)
# Box plot
#bp+scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest"))
# Scatter plot
#sp+scale_color_manual(values=wes_palette(n=3, name="GrandBudapest"))

labels <-j_fil_THP %>% 
  group_by(treat) %>% 
  summarise(N=paste0("n =", n()))
geom_text(data = labels,aes(x=treat,y=5e+5,label=N), size=6) 

labels_yH2AX_THP1_ <-High_spots_yH2AX_THP %>% 
  group_by(treat) %>% 
  summarise(N=paste0("n =", n()))

#geom_text(data = labels,aes(x=treat,y=5e+5,label=N), size=6) + 

#create biolin plot for yH2AX 
v_yH2AX <- ggplot(High_spots_yH2AX_THP, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+5)) +
  geom_text(data = labels_yH2AX_THP1_, aes(x=treat,y=5e+4,label=N), size=6) + 
  #stat_summary(fun.data = n_fun, geom = "text") +
  #scale_fill_few(name="Treatment") +
  ylab("yH2AX Nuclear Intensity Mean log10(A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("THP-1") + 
  labs(fill="Treatment") 
 

#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = High_spots_yH2AX_THP, paired = FALSE)
my_comparisons = list( c("30", "DMSO") )
stat_compare_means(comparisons = my_comparisons)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = High_spots_yH2AX_THP)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(High_spots_yH2AX_THP$nuc_intyH2AX_mean,    High_spots_yH2AX_THP$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")

#print violin plot with boxplot and significance values 
v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text= element_text(size=15))

scatter_THP <-ggplot(j_fil_THP, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) + 
  geom_point(shape=18, color="blue")+
  ggtitle("THP-1") + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")

scatter_THP<-ggplot(j_fil_THP, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=treat)) + 
  geom_point(alpha = 0.5)+
  ggtitle("THP-1") + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug() +
  scale_x_log10() + 
  scale_y_log10() + 
  geom_smooth(method= lm, se=FALSE, fullrange=TRUE)  

etp_THP <- filter(j_fil_THP, treat == "30")

etp_THP_Reg <- lm(nuc_intyH2AX_mean ~ spot_int53BP1_sum, data = etp_THP) 
summary(etp_THP_Reg)

summary(etp_THP_Reg)$r.squared



unt_THP <- filter(j_fil_THP, treat == "DMSO")

unt_THP_Reg <- lm(nuc_intyH2AX_mean ~ spot_int53BP1_sum, data = unt_THP) 
summary(unt_THP_Reg)

summary(unt_THP_Reg)$r.squared



scatter_THP + labs(fill="Treatment")
