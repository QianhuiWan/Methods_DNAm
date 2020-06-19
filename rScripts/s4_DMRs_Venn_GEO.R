
# load packages and data

library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(plyranges)
library(grid)
library(VennDiagram)
library(eulerr)
library(gridExtra)
# bumphunterDefault, ProbeLassoDefault, DMRcateDefault
load(file = here("Method_DNAme/rds/DMRmethods_DefaultPara_GEO_Grange.RData"))

# bumphunterOpti, ProbeLassoOpti, DMRcateOpti
load(file = here("Method_DNAme/rds/DMRmethods_OptiPara_GEO_Grange.RData"))

# bumphunterComp, ProbeLassoComp, DMRcateComp
load(file = here("Method_DNAme/rds/DMRmethods_ComparablePara_GEO_Grange.RData"))

##########################################################################################
# default para
## bumphunter, 99
## Probe Lasso, 2
## DMRcate, 395

find_overlaps(bumphunterDefault, ProbeLassoDefault) #0, 1_2
find_overlaps(bumphunterDefault, DMRcateDefault) #39, 1_3
find_overlaps(ProbeLassoDefault, DMRcateDefault) #2, 2_3

find_overlaps(find_overlaps(bumphunterDefault, ProbeLassoDefault), DMRcateDefault) #0

# Venn diagram
Venn_def <- draw.triple.venn(area1 = 99,
                             area2 = 2,
                             area3 = 395,
                             n12 = 0,
                             n13 = 39,
                             n23 = 2,
                             n123 = 0,
                             category = c("Bumphunter", "Probe Lasso", "DMRcate"),
                             fill = c("pink", "lightblue", "lightgoldenrodyellow"),
                             cex = 1.2, cat.cex = 0.9, 
                             cat.fontface = 2,
                             cat.just =list(c(0, -1.5), c(0.5, -3.5), c(1.5, -1.5)),
                             # lty = "blank"
                             cat.dist = c(0.05, 0.05, 0.05)
) 

Venn_def <- grid.arrange(gTree(children=Venn_def), vp=viewport(width = 0.9, height = 0.9),
                         top=textGrob("Default parameters",
                                      gp=gpar(fontface="bold",cex=0.8),
                                      x = 0.5, vjust = -0.8)) %>% 
  as_ggplot()

##############################################################################
# optimised para
## bumphunter, 2
## Probe Lasso, 133
## DMRcate, 691

find_overlaps(bumphunterOpti, ProbeLassoOpti) # 1, 1_2
find_overlaps(bumphunterOpti, DMRcateOpti) #1, 1_3
find_overlaps(ProbeLassoOpti, DMRcateOpti) #74, 2_3

find_overlaps(find_overlaps(bumphunterOpti, ProbeLassoOpti), DMRcateOpti) # 1

# Venn diagram
Venn_opti <- draw.triple.venn(area1 = 2,
                             area2 = 133,
                             area3 = 691,
                             n12 = 1,
                             n13 = 1,
                             n23 = 74,
                             n123 = 1,
                             category = c("Bumphunter", "Probe Lasso", "DMRcate"),
                             fill = c("pink", "lightblue", "lightgoldenrodyellow"),
                             cex = 1.2, cat.cex = 0.9, 
                             cat.fontface = 2,
                             cat.just =list(c(0, -0.5), c(1, -0.5), c(0.5, -1)),
                             # lty = "blank"
                             cat.dist = c(0.05, 0.05, 0.05)
) 

Venn_opti <- grid.arrange(gTree(children=Venn_opti), vp=viewport(width = 0.9, height = 0.9),
                         top=textGrob("Optimised parameters",
                                      gp = gpar(fontface="bold",cex=0.8),
                                      x = 0.5, vjust = -0.8)) %>% 
  as_ggplot()

# Comparable para ##################################################################
## bumphunter, 3
## Probe Lasso, 133
## DMRcate, 691

find_overlaps(bumphunterComp, ProbeLassoComp) #1, 1_2
find_overlaps(bumphunterComp, DMRcateComp) #1, 1_3
find_overlaps(ProbeLassoComp, DMRcateComp) #74 , 2_3

find_overlaps(find_overlaps(bumphunterComp, ProbeLassoComp), DMRcateComp) #0

Venn_comp <- draw.triple.venn(area1 = 3,
                              area2 = 133,
                              area3 = 691,
                              n12 = 1,
                              n13 = 1,
                              n23 = 74,
                              n123 = 1,
                              category = c("Bumphunter", "Probe Lasso", "DMRcate"),
                              fill = c("pink", "lightblue", "lightgoldenrodyellow"),
                              cex = 1.2, cat.cex = 0.9, 
                              cat.fontface = 2,
                              cat.just =list(c(0, -0.5), c(1, -0.5), c(0.5, -1)),
                              # lty = "blank"
                              cat.dist = c(0.05, 0.05, 0.05)) 

Venn_comp <- grid.arrange(gTree(children=Venn_comp), vp=viewport(width = 0.9, height = 0.9),
                          top=grid::textGrob("Comparable parameters",
                                        gp=gpar(fontface="bold",cex=0.8),
                                        x = 0.5, vjust = -0.8)) %>% 
  as_ggplot()

# jpeg(file = here("Method_DNAme/figures/test.jpeg"),
#      width = 8, height = 8, units = "cm", res = 320)
# Venn_comp
# dev.off()

VennPlots_fig <- ggarrange(Venn_def, Venn_opti, Venn_comp,
         labels = c("A", "B", "C", NULL),
         font.label = list(size = 12, face = "plain", color ="black"), 
         nrow=2, ncol = 2)

jpeg(file = here("Method_DNAme/figures/VennPlots_DMRs_v2.jpeg"),
     width = 16, height = 16, units = "cm", res = 320)
VennPlots_fig
dev.off()

