#The purpose of this script is to plot figure 2c 
#load libraries 
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
library(readxl)

#set working directory 
setwd("/Users/gallantkl")

#load datasets 
X_122_jurkat <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190218-121-KG-jurkat-yH2AX-53bp1-dose_20190218_141455[3607]/190221-EXP122-DoseResponse-Jurkat.xlsx")
X_119_BJAB <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190124-119-KLG-Jurkat-BJAB-THP1-Ad2x_20190124_113931[3510]/190218-119-BJAB.xlsx")
X_122_THP1 <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190218-122-KG-thp1-yH2AX-53bp1-dose_20190218_150512[3608]/190218-EXP122-DoseResponse-THP1.xlsx")

#prepare datasets 
#filter dataset to specific density 
j_fil_Jur <- X_122_jurkat %>%
  filter(treat %in% c("DMSO", "30")) 

j_fil_BJAB <- X_119_BJAB %>%
  filter(treat %in% c("DMSO", "30"))

j_fil_THP <- X_122_THP1 %>%
  filter(treat %in% c("DMSO", "30"))


#filter high yh2AX 
high_yH2AX_Jur <- filter(j_fil_Jur, treat_yH2AX == "High")
high_spots_yH2AX2 <- filter(j_fil_BJAB, treat_yH2AX == "High")
High_spots_yH2AX_THP <- filter(j_fil_THP, treat_yH2AX == "High")
#plot 

#create my comparisons 
my_comparisons = list( c("30", "DMSO") )

#plot jurkat ----------
labels <-j_fil_Jur %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))

v_yH2AX <- ggplot(j_fil_Jur, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels,aes(x=treat,y=1e+5,label=N), size=6) + 
  #stat_summary(fun.data = n_fun, geom = "text") +
  #scale_fill_few(name="Treatment") +
  #labs(fill = "Treatment") + 
  ylab(expression(gamma ~ "H2AX Nuclear Intensity Mean (A.U) + 1")) + 
  #xlab("Treatment") + 
  ggtitle("Jurkat")  


#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = High_spots_Jurkat, paired = FALSE)
stat_compare_means(comparisons = my_comparisons)

#print violin plot with boxplot and significance values 
fig2c_jurkat <- v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  
  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + 
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())


#plot bjab ---------------

labels <-j_fil_BJAB %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))

v_yH2AX <- ggplot(j_fil_BJAB, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels,aes(x=treat,y=1e+5,label=N), size=6) + 
  #stat_summary(fun.data = n_fun, geom = "text") +
  #scale_fill_few(name="Treatment") +
  #labs(fill = "Treatment") + 
  ylab(expression(gamma ~ "H2AX Nuclear Intensity Mean (A.U) + 1")) + 
  #xlab("Treatment") + 
  ggtitle("BJAB")  


#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = high_spots_yH2AX2, paired = FALSE)
stat_compare_means(comparisons = my_comparisons)

#print violin plot with boxplot and significance values 
fig2c_bjab <- v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  
  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + 
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())


#plot thp1 ------------

labels <-j_fil_THP %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))

v_yH2AX <- ggplot(j_fil_THP, aes(x= treat, y= nuc_intyH2AX_mean +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+1, 1e+6)) +
  geom_text(data = labels,aes(x=treat,y=1e+5,label=N), size=6) + 
  #stat_summary(fun.data = n_fun, geom = "text") +
  #scale_fill_few(name="Treatment") +
  #labs(fill = "Treatment") + 
  ylab(expression(gamma ~ "H2AX Nuclear Intensity Mean (A.U) + 1")) + 
  #xlab("Treatment") + 
  ggtitle("THP-1")  


#determine significance values for v_yH2AX 
compare_means(nuc_intyH2AX_mean ~ treat,  data = High_spots_yH2AX_THP, paired = FALSE)
stat_compare_means(comparisons = my_comparisons)

#print violin plot with boxplot and significance values 
fig2c_thp1 <- v_yH2AX + geom_boxplot(width=0.1, notch = TRUE, show.legend = FALSE) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  
  scale_fill_brewer(palette = "BuGn") + theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + 
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())



#--arrange plots--
fig2c <- ggarrange(fig2c_jurkat, fig2c_bjab, fig2c_thp1,
                    ncol = 3, nrow = 1)
fig2c
