#the purpose of this code is for data visualization of cells 
library(tidyverse)
library(ggthemes)
library(scales)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(stats)
library(forcats)
library(readxl)

#set working directory 
setwd("/Users/gallantkl")

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


#filter high 53bp1
High_spots_Jurkat <- filter(j_fil_Jur, treat_lev == "Positive")
High_spots_BJAB <- filter(j_fil_BJAB, treat_lev == "Positive")
High_spots_THP <- filter(j_fil_THP, treat_lev == "High")
#plot 

#create my comparisons 
my_comparisons = list( c("30", "DMSO") )
#filter dataset to high spots 


#plot jurkat -------
labels <-j_fil_Jur %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))

v_53BP1 <- ggplot(j_fil_Jur, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+9)) +
  geom_text(data = labels,aes(x=treat,y=1e+8,label=N), size=6) + 
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  ggtitle("Jurkat") + 
  annotate("text", x=1.35, y=500000, label = "33.4%", size = 6) +  #untreated and high
  annotate("text", x=1.35, y=10, label = "66.7%", size = 6) + #untreated and low 
  annotate("text", x=2.45, y=500000, label = "69.0%", size = 6) + #ETP and high 
  annotate("text", x=2.45, y=10, label = "31.0%", size = 6)  #EtP and low 


#high spots stats
compare_means(spot_int53BP1_sum ~ treat,  data = High_spots_Jurkat, paired = FALSE)
my_comparisons = list( c("30", "DMSO") )
stat_compare_means(comparisons = my_comparisons)


#print violin plot 
fig2b_jurkat <- v_53BP1 + scale_fill_brewer(palette = "BuGn") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + 
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())

#--bjab plot 

labels <-j_fil_BJAB %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))

v_53BP1 <- ggplot(j_fil_BJAB, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+9)) +
  geom_text(data = labels,aes(x=treat,y=1e+8,label=N), size=6) + 
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  ggtitle("BJAB") + 
  annotate("text", x=1.35, y=50000, label = "5.5%", size = 6) +  #untreated and high
  annotate("text", x=1.35, y=10, label = "94.5%", size = 6) + #untreated and low 
  annotate("text", x=2.35, y=50000, label = "41.8%", size = 6) + #ETP and high 
  annotate("text", x=2.35, y=10, label = "58.2%", size = 6)  #EtP and low 



#high spots stats
compare_means(spot_int53BP1_sum ~ treat,  data = High_spots_BJAB, paired = FALSE)
my_comparisons = list( c("30", "DMSO") )
stat_compare_means(comparisons = my_comparisons)


#print violin plot 
fig2b_bjab <- v_53BP1 + scale_fill_brewer(palette = "BuGn") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + 
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())



#--thp1 plot 
labels <-j_fil_THP %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))

v_53BP1 <- ggplot(j_fil_THP, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) +
  geom_violin(show.legend = FALSE) + 
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+9)) +
  geom_text(data = labels,aes(x=treat,y=1e+8,label=N), size=6) + 
  scale_fill_few(name="Treatment") +
  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  ggtitle("THP-1") + 
  annotate("text", x=1.39, y=80000, label = "61.8%", size = 6) +  #untreated and high
  annotate("text", x=1.39, y=10, label = "38.2%", size = 6) + #untreated and low 
  annotate("text", x=2.42, y=80000, label = "71.2%", size = 6) + #ETP and high 
  annotate("text", x=2.42, y=10, label = "28.8%", size = 6)


#high spots stats
compare_means(spot_int53BP1_sum ~ treat,  data = High_spots_THP, paired = FALSE)
my_comparisons = list( c("30", "DMSO") )
stat_compare_means(comparisons = my_comparisons)


#print violin plot 
fig2b_thp1 <- v_53BP1 + scale_fill_brewer(palette = "BuGn") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  theme_bw(base_size = 20) + 
  scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) + 
  theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())


#---arrange
fig2b <- ggarrange(fig2b_jurkat, fig2b_bjab, fig2b_thp1,
                   ncol = 3, nrow = 1)
fig2b

