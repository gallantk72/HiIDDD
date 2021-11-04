#the purpose of this code is for data visualization of cells 
library(tidyverse)
library(ggthemes)
library(scales)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(stats)
library(forcats)
library(ggsignif)


#find and read data from the directory
library(readxl)
setwd("/Users/gallantkl")

X_122_jurkat <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190218-121-KG-jurkat-yH2AX-53bp1-dose_20190218_141455[3607]/190221-EXP122-DoseResponse-Jurkat.xlsx")
View(X_122_jurkat)


lev1 <- c("DMSO", "30")



#filter 0 and 30 etp and reorder 
j_fil_Jur <- X_122_jurkat %>%
  filter(treat %in% c("DMSO", "30")) 
  #mutate(names = factor(treat, levels = c("DMSO","30"))) %>% 

j_fil_Jur$treat
scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP"))
#reorder groups
#reorder groups 
j_fil_Jur$treat <- factor(j_fil_Jur$treat, levels = c("DMSO","30"))
#j_fil_Jur$treat <- as.factor(j_fil_Jur$treat)

high_yH2AX_Jur <- filter(j_fil_Jur, treat_yH2AX == "High")

table(high_yH2AX_Jur$treat)

#-------------------------------------------
#CHI SQUARE FOR yH2AX
#percentages 
j_Jur_perc <- j_fil_Jur %>%
  group_by(treat,treat_yH2AX ) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts))

#CHI SQUARE ANALYSIS 
obsfreq <- c(2317,328)
nullprobs <- c(0.5,0.5)
#chisq.test(obsfreq,p=nullprobs)

obsfreq2 <- c(954,700)
nullprobs <- c(0.5,0.5)

rows = 2 
matriz_Jur = matrix(c(obsfreq,obsfreq2), 
                    nrow=rows,
                    byrow=TRUE)

rownames(matriz_Jur) = c("DMSO","ETP")
colnames(matriz_Jur) = c("Low", "High")

matriz_Jur
chi_matriz_Jur <-chisq.test(matriz_Jur)
#________________________________

#CHI SQUARE FOR 53BP1 
j_fil_Jurperc <- j_fil_Jur %>%
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
matriz_Jur = matrix(c(obsfreq,obsfreq2), 
                     nrow=rows,
                     byrow=TRUE)

rownames(matriz_Jur) = c("DMSO","ETP")
colnames(matriz_Jur) = c("Low", "High")

matriz_Jur
chi_matriz_Jur <-chisq.test(matriz_Jur)

#--------------------
#filled contour plots 
scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP"))

Jurkat_confilled <- ggplot(j_fil_Jur, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat), show.legend = FALSE) +
  xlab("53BP1 Sum Spot Intensity (A.U.)") + 
  ylab(expression(gamma ~ "H2AX Mean Nuclear Intensity (A.U.)")) + 
  labs(fill="Treatment") + 
  ggtitle("Jurkat") + 
  scale_x_continuous(trans = 'log2', labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific)

#Jurkat_confilled +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
Jurkat_confilled + theme_bw(base_size = 20) +
  theme(axis.text = element_text(size=15)) +#, legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP"))



#Jur_confilled +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
Jur_confilled + theme_bw(base_size = 20)
#prep data for barplot
j_fil_Jurperc <- j_fil_Jur %>%
  group_by(treat, treat_lev) %>%
  summarise(counts = n()) %>%
  mutate(perc=counts/sum(counts)) 
  #ggplot( aes(x = treat, y = perc*100)) +
  #geom_bar(
    #aes(color = treat_lev, fill = treat_lev),
    #stat = "identity", position = position_stack()
  #) 
  #labs(x = "Treatment", y = "Percentage of Total Cells", fill = "Treatment Level") 
  
 
#mutate(name = fct_relevel(treat,
#"DMSO", "30")) %>% 




j_barplot <- ggplot(j_fil_Jurperc, aes(x = treat, y = perc*100)) +
  geom_bar(
    aes(fill = treat_lev),
    stat = "identity", position = position_stack()
  ) + 
  ggtitle("Jurkat") + 
  labs(x = "Treatment", y = "Total Cells (%)", fill = "53BP1 Spot Status", main = "Jurkat")

#print barplot
j_barplot + scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text = element_text(size = 15))

#ggplot(j_fil_Jur, aes(x = treat, y = perc*100)) +
 # geom_bar(aes(fill = treat_lev, stat = "identity")) +
 # labs(x = "Treatment", y = "Percentage of 53BP1 Positive Cells", fill = "53BP1 Positive")

geom_text(
  aes(label = perc*100, group = treat_lev)
)

names <- c("DMSO", "30")

lev1 <- c("DMSO", "30")

j_fil_cat <- j_fil %>%
  mutate(names = factor(names, levels = lev1))
  

#msleep %>% 
  #select(order, name, sleep_total) %>% 
  #filter(order %in% c("Didelphimorphia", "Diprotodontia"))
  #df <- tibble(names, values) %>%
  #mutate(names = factor(names, levels = lev1))

perc_pos <- filter(j_fil_Jurperc, treat_lev == "Positive")

#create function to call sample sizes 
n_fun <- function(x){
  return(data.frame(y = max(x)*1.1, 
                    label = paste0("n = ", perc_pos$perc)))
}

#stats 
High_spots_Jurkat <- filter(j_fil_Jur, treat_lev == "Positive")

High_Jur <- High_spots_Jurkat %>%
  group_by(treat) %>% 
  get_summary_stats(spot_int53BP1_sum, type="median_mad")

yH2AX_Jur <- j_fil_Jur %>%
  group_by(treat) %>% 
  get_summary_stats(nuc_intyH2AX_mean, type="median_mad")
#create violin plot for 53BP1 
v_53BP1 <- ggplot(High_spots_Jurkat, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(size = 0.75) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+8)) +
  stat_summary(fun.data = n_fun, geom = "text") +
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("Jurkat") + 
  labs(fill="Treatment") + 
  theme_bw()
#stats
  compare_means(spot_int53BP1_sum ~ treat,  data = High_spots_Jurkat, paired = FALSE)
  my_comparisons = list( c("30", "DMSO") )
  stat_compare_means(comparisons = my_comparisons)
  
  #perform one way ANOVA on yH2AX data and print summary 
  ad_aov <- aov(spot_int53BP1_sum ~ treat, data = High_spots_Jurkat)
  summary(ad_aov)
  
  #post-hoc
  ad_pairwise <- pairwise.t.test(High_spots_Jurkat$spot_int53BP1_sum,    High_spots_Jurkat$treat, p.adj = "none")
  # tidy the post hoc
  tidy_ad_pairwise <- broom::tidy(ad_pairwise)
  
  # look at the comparisons and p-values
  head(tidy_ad_pairwise)
  #print violin plot with stats 
  v_53BP1 + geom_boxplot(width=0.1, notch = TRUE, boxlwd = 2) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") + scale_fill_brewer(palette = "BuGn") + theme_bw() + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP"))
  
  
#print violin plot 
v_53BP1 + scale_fill_brewer(palette = "BuGn") + theme_bw()

#j_fil_Jur$treat <- factor(j_fil_Jur$treat, levels = c("DMSO","ETP"))
#create function to call sample sizes 
n_fun <- function(x){
  return(data.frame(y = max(x)*1.1, 
                   label = perc*100, data=perc_pos))
}

#n.fun
#n_fun <- function(x){
  #return(data.frame(y= quantile(y,probs=0.95)*1.1, 
                    #label = paste0("n = ",length(x)), size = 6))
#}

labels <-j_fil_Jur %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))

#geom_text(data = labels,aes(x=treat,y=5e+7,label=N), size=6) + 

#X position
#Labels$mpg <- 25
#Plot
#ggplot(mtcars, aes(x = mpg, stat = "count",fill=as.factor(carb))) + geom_histogram(bins = 8)+

#paste0("n = ", perc_pos$perc
#create violin plot for 53BP1 positive 53BP1
v_53BP1 <- ggplot(j_fil_Jur, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+9)) +
  geom_text(data = labels,aes(x=treat,y=1e+8,label=N), size=6) + 
  #geom_text(aes(label=..count..), y=1000000, stat='count', size=6) + 
  #stat_summary(fun.data = n_fun, geom = "text", hjust =0.5, vjust=0.9) +
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum log10(A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("Jurkat") + 
  labs(fill="Treatment") + 
  annotate("text", x=1.35, y=500000, label = "33.4%", size = 6) +  #untreated and high
  annotate("text", x=1.35, y=10, label = "66.7%", size = 6) + #untreated and low 
  annotate("text", x=2.45, y=500000, label = "69.0%", size = 6 ) + #ETP and high 
  annotate("text", x=2.45, y=10, label = "31.0%", size = 6)  #EtP and low 
#theme_bw()

#high-spots 
High_spots_Jurkat <- filter(j_fil_Jur, treat_lev == "Positive")
#High_spots_Jurkat$treat <- factor(High_spots_Jurkat$treat, levels = c("DMSO", "ETP"))
#high spots stats
compare_means(spot_int53BP1_sum ~ treat,  data = High_spots_Jurkat, paired = FALSE)
my_comparisons = list( c("30", "DMSO") )
stat_compare_means(comparisons = my_comparisons)

#print summary data
aggregate(spot_int53BP1_sum ~ treat, High_spots_Jurkat,mean)
aggregate(spot_int53BP1_sum ~ treat, High_spots_Jurkat,sd)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(spot_int53BP1_sum ~ treat, data = High_spots_Jurkat)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(High_spots_Jurkat$spot_int53BP1_sum,    High_spots_Jurkat$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#reorder
#j_fil_Jur$treat <- factor(j_fil_Jur$treat, levels = c("DMSO","ETP"))

#print violin plot 
v_53BP1 + scale_fill_brewer(palette = "BuGn") + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text = element_text(size = 15))
                                                                                                                                         
                                                                                                                                                      

# Box plot
#bp+scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest"))
# Scatter plot
#sp+scale_color_manual(values=wes_palette(n=3, name="GrandBudapest"))

#create biolin plot for yH2AX 
v_yH2AX <- ggplot(j_fil_Jur, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels,aes(x=treat,y=1e+5,label=N), size=6) + 
  #stat_summary(fun.data = n_fun, geom = "text") +
  scale_fill_few(name="Treatment") +
  labs(fill = "Treatment") + 
  ylab(expression(gamma ~ "H2AX Nuclear Intensity Mean log10(A.U) + 1")) + 
  #xlab("Treatment") + 
  ggtitle("Jurkat")  


#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = j_fil_Jur, paired = FALSE)
my_comparisons = list( c("30", "DMSO") )
stat_compare_means(comparisons = my_comparisons)

#summary
aggregate(spot_int53BP1_sum ~ treat, High_spots_Jurkat,mean)
aggregate(spot_int53BP1_sum ~ treat, High_spots_Jurkat,sd)
aggregate(spot_int53BP1_sum ~ treat, j_fil_Jur,mean)
aggregate(spot_int53BP1_sum ~ treat, j_fil_Jur,sd)
aggregate(nuc_intyH2AX_mean ~ treat, j_fil_Jur,mean)
aggregate(nuc_intyH2AX_mean ~ treat, j_fil_Jur,sd)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = j_fil_Jur)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(j_fil_Jur$nuc_intyH2AX_mean,    j_fil_Jur$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")
#geom_signif(comparisons = list(c("G1","G2")), y_position = 28,
           # tip_length = 0, vjust = .1
#print violin plot with boxplot and significance values 
v_yH2AX + 
  geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  
  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + 
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element())


#v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, y_position = 1e+4) +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP"))

scatter_jur <-ggplot(j_fil_Jur, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) + 
  geom_point(shape=18, color="blue")+
  ggtitle("Jurkat") + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")

scatter_jur<-ggplot(j_fil_Jur, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean, color=treat)) + 
  geom_point(alpha = 0.5)+
  ggtitle("Jurkat") + 
  xlab("53BP1 Sum Spot Intensity") + 
  ylab("yH2AX Mean Nuclear Intensity") + 
  geom_rug() +
  scale_x_log10() + 
  scale_y_log10() + 
  geom_smooth(method= lm, se=FALSE, fullrange=TRUE)  

etp_Jur <- filter(j_fil_Jur, treat == "30")

etp_Jur_Reg <- lm(nuc_intyH2AX_mean ~ spot_int53BP1_sum, data = etp_Jur) 
summary(etp_Jur_Reg)

summary(etp_Jur_Reg)$r.squared



unt_Jur <- filter(j_fil_Jur, treat == "DMSO")

unt_Jur_Reg <- lm(nuc_intyH2AX_mean ~ spot_int53BP1_sum, data = unt_Jur) 
summary(unt_Jur_Reg)

summary(unt_Jur_Reg)$r.squared



#scatter_jur + facet_grid(vars(treat)) + labs(fill="Treatment") 

scatter_jur + labs(fill="Treatment") 

View(high_yH2AX_Jur)
#---------

labels_yH2AX_Jur <- high_yH2AX_Jur %>% 
  group_by(treat) %>% 
  summarise(N=paste0("n =", n()))

v_yH2AX <- ggplot(high_yH2AX_Jur, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+5)) +
  geom_text(data = labels_yH2AX_Jur,aes(x=treat,y=5e+4,label=N), size=6) + 
  scale_fill_few(name="Treatment") +
  labs(fill = "Treatment") + 
  ylab("Pos-yH2AX Nuclear Intensity Mean log10(A.U) + 1") + 
  xlab("Treatment") + 
  ggtitle("Jurkat")  


#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = high_yH2AX_Jur, paired = FALSE)
my_comparisons = list( c("30", "DMSO") )
stat_compare_means(comparisons = my_comparisons)

#summary
aggregate(spot_int53BP1_sum ~ treat, High_spots_Jurkat,mean)
aggregate(spot_int53BP1_sum ~ treat, High_spots_Jurkat,sd)
aggregate(spot_int53BP1_sum ~ treat, j_fil_Jur,mean)
aggregate(spot_int53BP1_sum ~ treat, j_fil_Jur,sd)
aggregate(nuc_intyH2AX_mean ~ treat, j_fil_Jur,mean)
aggregate(nuc_intyH2AX_mean ~ treat, j_fil_Jur,sd)

#perform one way ANOVA on yH2AX data and print summary 
ad_aov <- aov(nuc_intyH2AX_mean ~ treat, data = j_fil_Jur)
summary(ad_aov)

#post-hoc
ad_pairwise <- pairwise.t.test(j_fil_Jur$nuc_intyH2AX_mean,    j_fil_Jur$treat, p.adj = "none")
# tidy the post hoc
tidy_ad_pairwise <- broom::tidy(ad_pairwise)

# look at the comparisons and p-values
head(tidy_ad_pairwise)

#save to excel
#write.csv(tidy_ad_pairwise, "tidy_ad_pairwise_20601.csv")
#geom_signif(comparisons = list(c("G1","G2")), y_position = 28,
# tip_length = 0, vjust = .1
#print violin plot with boxplot and significance values 
v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + theme(axis.text = element_text(size = 15))


