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
library(rstatix)
library(PMCMRplus)
library(rstatix)
library(readxl)
library(plotly)
#setwd("/Users/gallantkl")

#exp206_all <- read_excel("Desktop/exp206_CD4_KG/exp206_allrep_CD4_KG.xlsx")

#exp206_all <- read_excel("data/raw-data_v1/fig_4/exp206_allrep_CD4_KG.xlsx")

exp206_all <- list.files(path = "data/raw-data_v1/fig_4/exp206-files",    
                       pattern = "*.xlsx", full.names = TRUE) %>% 
  lapply(read_excel) %>%                                           
  bind_rows                                                      

exp206_all
View(exp206_all)

#create dataframe for means (53BP1)
repli_data <- exp206_all %>%
  group_by(treat)
# get_summary_stats(spot_int53BP1_sum, type="mean_sd")
#View(repli_data)
#reorder 
repli_data <- exp206_all %>%
  group_by(treat, replicate) %>% 
  summarise_all(mean) 

repli_data$treat <- factor(repli_data$treat, levels = c("Untreated","DMSO","ETP"))

my.comparisons = list( c("Untreated", "DMSO"), c("DMSO","ETP"), c("ETP","Untreated") ) 
#create dotplot for 53BP1
repli_53BP1 <- ggplot(repli_data, aes(x = treat, y = spot_int53BP1_sum)) +
  geom_dotplot(aes(fill = treat), binaxis = "y", binwidth = 0.1, stackdir= "center", show.legend = FALSE, dotsize = 1.4) + 
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.1) +
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+7)) + 
  ylab("53BP1 spot intensity sum (a.u.) + 1") + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("CD4+T cells") + 
  scale_fill_brewer(palette = "BuGn") + 
  theme_classic(base_size = 20) + 
  theme(axis.text = element_text(size = 15)) + 
  #theme_bw(base_size = 20) + 
  stat_compare_means(comparisons = my.comparisons, label = "p.signif", method = "t.test", paired = FALSE, ) +    
  annotate("text", x=2.8, y=2, label = "ANOVA, p = 0.99", size = 6) 
  #theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())

repli_53BP1
#dataframe for means yH2AX 

#create comparisons for t-test
my.comparisons.2 = list( c("Untreated", "DMSO"), c("DMSO","ETP"), c("ETP","Untreated") )


#ggplot dotplot 
repli_H2AX <- ggplot(repli_data, aes(x = treat, y = nuc_intyH2AX_mean)) +
  geom_dotplot(aes(fill = treat), binaxis = "y", binwidth = 0.1, stackdir= "center", show.legend = FALSE, dotsize = 1.4) + 
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.1) +
  scale_y_continuous(trans = "log10", limits = c(1e+0, 1e+7)) + 
  ylab(expression(gamma ~ "H2AX nuclear intensity mean (a.u) + 1")) + 
  xlab("Treatment") + 
  labs(fill="Treatment") + 
  ggtitle("CD4+T cells") + 
  scale_fill_brewer(palette = "BuGn") + 
  theme_classic(base_size = 20) + 
  #theme_bw(base_size = 20) + 
  #stat_pvalue_manual(spws,label ="p.adj.signif") +  
  stat_compare_means(comparisons = my.comparisons.2, label = "p.signif", method = "t.test", paired = FALSE) +      
  annotate("text", x=2.8, y=2, label = "ANOVA, p = 0.93", size = 6 ) 
  #theme(axis.text = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())  #t.test comparison 

#ANOVA 
compare_means(mean ~ repli, data = repli_data_H2AX, method = "anova") 

my_comparisons = list( c("Untreated", "DMSO"), c("DMSO","ETP"), c("ETP","Untreated") ) 
stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") 

repli_H2AX  #+ stat_compare_means(comparisons = my_comparisons, label = "p.signif", paired= FALSE)

#create figure 
fig4a <- ggarrange(repli_53BP1,repli_H2AX, 
                   ncol = 2, nrow = 1)
fig4a
