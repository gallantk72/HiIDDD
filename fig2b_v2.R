# the purpose of this R script is to load all graphs from figures 2 

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
library(readxl)

#set working directory 
setwd("/Users/gallantkl")

#load datasets 
X_122_jurkat <- read_excel("Desktop/190218-121-122-THP1-Jurkat-ETPtitration-BJAB/190218-121-KG-jurkat-yH2AX-53bp1-dose_20190218_141455[3607]/190221-EXP122-DoseResponse-Jurkat.xlsx")


#prepare datasets 
#filter dataset to specific density 
j_fil_Jur <- X_122_jurkat %>%
  filter(treat %in% c("DMSO", "30")) 

#filter dataset to high spots 
High_spots_Jurkat <- filter(j_fil_Jur, treat_lev == "Positive")

#create violin plots for fig 2b (53BP1)

#create 53BP1 format function 
violin_53BP1 <-  ylab("53BP1 Spot Intensity Sum (A.U) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
  theme_bw(base_size = 20) + 
  theme(axis.text = element_text(size = 15))


#define size of each group 
labels <-j_fil_Jur %>% 
  group_by(treat) %>% 
  dplyr::summarise(N=paste0("n =", n()))

#create custom annotations for 53BP1 plots 
jurkat_ann <- annotate("text", x=1.35, y=500000, label = "33.4%", size = 6) +  #untreated and high
  annotate("text", x=1.35, y=10, label = "66.7%", size = 6) + #untreated and low 
  annotate("text", x=2.45, y=500000, label = "69.0%", size = 6 ) + #ETP and high 
  annotate("text", x=2.45, y=10, label = "31.0%", size = 6)  #ETP and low 


#define my_comparisons for stats test   
my_comparisons <- list( c("30", "DMSO") )

#define aes for each cell type 
jurkat_53BP1 <- ggplot(j_fil_Jur, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) + ggtitle("Jurkat") + geom_text(data = labels,aes(x=treat,y=1e+8,label=N), size=6) + scale_fill_few(name="Treatment") +
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+9)) + scale_fill_brewer(palette = "BuGn") + scale_x_discrete(labels = c("DMSO"= "DMSO", "30"= "ETP")) 
bjab_53BP1 <- ggplot(j_fil_BJAB, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat)) 
thp1_53BP1 <- ggplot(j_fil_THP, aes(x= treat, y= spot_int53BP1_sum +1, fill= treat))


#Plot individual plots 
compare_means(spot_int53BP1_sum ~ treat,  data = High_spots_Jurkat, paired = FALSE)

jurkat_53BP1 + violin_53BP1 + jurkat_ann 








mybarplot <- function(mydf, myxcol, myycol, mytitle) {
  ggplot2::ggplot(data = mydf, aes(x=reorder({{ myxcol }}, 
                                             {{ myycol }}), y= {{ myycol }})) +
    geom_col(color = "black", fill="#0072B2") +
    xlab("") +
    ylab("") +
    coord_flip() +
    ggtitle(mytitle) +
    theme_classic()   +
    theme(plot.title=element_text(size=24))
}

Jurkat_confilled <- ggplot(j_fil_Jur, aes(x= spot_int53BP1_sum, y= nuc_intyH2AX_mean)) + stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = treat), show.legend = FALSE) 
 
  
confilled_format <- xlab("53BP1 Sum Spot Intensity (A.U.)") + 
 ylab(expression(gamma ~ "H2AX Mean Nuclear Intensity (A.U.)")) + 
  labs(fill="Treatment") + 
  ggtitle("Jurkat") + 
  scale_x_continuous(trans = 'log2', labels =scales::scientific) + 
  scale_y_continuous(trans = 'log2', labels =scales::scientific) +
  theme_bw(base_size = 20) + 
  theme(axis.text = element_text(size=15), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

  

#Jurkat_confilled +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

mybarplot <- function(mydf, myxcol, myycol, mytitle) {
  ggplot2::ggplot(data = mydf, aes(x=reorder({{ myxcol }}, 
                                             {{ myycol }}), y= {{ myycol }})) +
    geom_col(color = "black", fill="#0072B2") +
    xlab("") +
    ylab("") +
    coord_flip() +
    ggtitle(mytitle) +
    theme_classic()   +
    theme(plot.title=element_text(size=24))
  
  
}

