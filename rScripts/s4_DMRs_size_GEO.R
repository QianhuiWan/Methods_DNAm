
library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(plyranges)
library(minfi)

# bumphunterDefault, ProbeLassoDefault, DMRcateDefault
load(file = here("Method_DNAme/rds/DMRmethods_DefaultPara_GEO_Grange.RData"))

# bumphunterOpti, ProbeLassoOpti, DMRcateOpti
load(file = here("Method_DNAme/rds/DMRmethods_OptiPara_GEO_Grange.RData"))

# bumphunterComp, ProbeLassoComp, DMRcateComp
load(file = here("Method_DNAme/rds/DMRmethods_ComparablePara_GEO_Grange.RData"))


# default para ###########################################
## bumphunter
bumphunter_def <- bumphunterDefault %>% width() %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("Bumphunter", nrow(.))) %>% 
  magrittr::inset("Parameter", value=rep("Default", nrow(.)))

## Probe Lasso
ProbeLasso_def <- ProbeLassoDefault %>% width() %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("Probe Lasso", nrow(.))) %>% 
  magrittr::inset("Parameter", value=rep("Default", nrow(.)))

## DMRcate
DMRcate_def <- DMRcateDefault %>% width() %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("DMRcate", nrow(.))) %>% 
  magrittr::inset("Parameter", value=rep("Default", nrow(.)))

# optimised para ############################################
## bumphunter
bumphunter_opti <- bumphunterOpti %>% width() %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("Bumphunter", nrow(.))) %>% 
  magrittr::inset("Parameter", value=rep("Optimised", nrow(.)))

## Probe Lasso
ProbeLasso_opti <- ProbeLassoOpti %>% width() %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("Probe Lasso", nrow(.))) %>% 
  magrittr::inset("Parameter", value=rep("Optimised", nrow(.)))

## DMRcate
DMRcate_opti <- DMRcateOpti %>% width() %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("DMRcate", nrow(.))) %>% 
  magrittr::inset("Parameter", value=rep("Optimised", nrow(.)))


# Comparable para ##################################################
## bumphunter
bumphunter_comp <- bumphunterComp %>% width() %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("Bumphunter", nrow(.))) %>% 
  magrittr::inset("Parameter", value=rep("Comparable", nrow(.)))

## Probe Lasso
ProbeLasso_comp <- ProbeLassoComp %>% width() %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("Probe Lasso", nrow(.))) %>% 
  magrittr::inset("Parameter", value=rep("Comparable", nrow(.)))

## DMRcate
DMRcate_comp <- DMRcateComp %>% width() %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("DMRcate", nrow(.))) %>% 
  magrittr::inset("Parameter", value=rep("Comparable", nrow(.)))



allWidth <- rbind(bumphunter_def, ProbeLasso_def, DMRcate_def, 
                  bumphunter_opti, ProbeLasso_opti, DMRcate_opti,
                  bumphunter_comp, ProbeLasso_comp, DMRcate_comp) %>% 
  `colnames<-`(c("DMR length", "Methods", "Parameters"))

# rbind all rows together
theme_set(theme_pubr())
DMR_size_plot <- allWidth %>% 
  ggplot(aes(x=factor(Methods, levels = c("Bumphunter", "Probe Lasso", "DMRcate")), 
             y=`DMR length`, 
             fill=factor(Parameters, levels = c("Default", "Optimised", "Comparable"))))+
  geom_boxplot()+
  # facet_grid(.~Methods)+
  # geom_violin()+
  # stat_compare_means(label = "p.signif", method = "wilcox.test", size = 5)+
  # geom_signif(comparisons = list(1:3), map_signif_level=TRUE)+ #,textsize = 6 , test = "wilcox.test"
  labs(y="DMR length", x="", fill="Parameters")+
  theme(
    text = element_text(size=12, colour = "black"),
    axis.text.x = element_text(size=12, colour = "black"),
    axis.title.x = element_text(size=12, colour = "black"),
    axis.text.y=element_text(size=12, colour = "black"), 
    legend.position = "right")

ggsave(filename = "DMR_size_default&Opti_v3.jpeg",
       plot = DMR_size_plot,
       device = "jpeg",
       path = here("Method_DNAme/rds_figure/"),scale = 1.6 ,
       width = 17, height = 6, units = "cm", dpi = 320)


