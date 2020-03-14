
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
load(file = here("Method_DNAme/rds/DMRmethods_DefaultPara_Grange.RData"))

# bumphunterOpti, ProbeLassoOpti, DMRcateOpti
load(file = here("Method_DNAme/rds/DMRmethods_OptiPara_Grange.RData"))

# bumphunterComp, ProbeLassoComp, DMRcateComp
load(file = here("Method_DNAme/rds/DMRmethods_ComparablePara_Grange.RData"))

##############################################################################
# default para
## bumphunter, 23
## Probe Lasso, 0
## DMRcate, 1199

find_overlaps(bumphunterDefault, ProbeLassoDefault) #0, 1_2
find_overlaps(bumphunterDefault, DMRcateDefault) #21, 1_3
find_overlaps(ProbeLassoDefault, DMRcateDefault) #0, 2_3

find_overlaps(find_overlaps(bumphunterDefault, ProbeLassoDefault), DMRcateDefault) #0

# Venn diagram
Venn_def <- draw.triple.venn(area1 = 23,
                             area2 = 0,
                             area3 = 1199,
                             n12 = 0,
                             n13 = 21,
                             n23 = 0,
                             n123 = 0,
                             category = c("Bumphunter", "Probe Lasso", "DMRcate"),
                             fill = c("pink", "lightblue", "lightgoldenrodyellow"),
                             cex = 1.2, cat.cex = 0.6, 
                             cat.fontface = 2,
                             cat.just =list(c(0, -0.5), c(1, -0.5), c(0.5, -1.5)),
                             # lty = "blank"
                             cat.dist = c(0.05, 0.05, 0.05)
) 

Venn_def <- grid.arrange(gTree(children=Venn_def), 
                         top=textGrob("Default parameters",
                                      gp=gpar(fontface="bold",cex=0.8))) %>% 
  as_ggplot()


##############################################################################
# optimised para
## bumphunter, 0
## Probe Lasso, 122
## DMRcate, 1045

find_overlaps(bumphunterOpti, ProbeLassoOpti) # 0, 1_2
find_overlaps(bumphunterOpti, DMRcateOpti) #0, 1_3
find_overlaps(ProbeLassoOpti, DMRcateOpti) #67, 2_3

find_overlaps(find_overlaps(bumphunterOpti, ProbeLassoOpti), DMRcateOpti) # 0

# Venn diagram
Venn_opti <- draw.triple.venn(area1 = 0,
                             area2 = 122,
                             area3 = 1045,
                             n12 = 0,
                             n13 = 0,
                             n23 = 67,
                             n123 = 0,
                             category = c("Bumphunter", "Probe Lasso", "DMRcate"),
                             fill = c("pink", "lightblue", "lightgoldenrodyellow"),
                             cex = 1.2, cat.cex = 0.6, 
                             cat.fontface = 2,
                             cat.just =list(c(0, -0.5), c(1, -0.5), c(0.5, -1.5)),
                             # lty = "blank"
                             cat.dist = c(0.05, 0.05, 0.05)
) 

Venn_opti <- grid.arrange(gTree(children=Venn_opti), 
                         top=textGrob("Optimised parameters",
                                      gp=gpar(fontface="bold",cex=0.8))) %>% 
  as_ggplot()


# Comparable para ##################################################################
## bumphunter, 2
## Probe Lasso, 122
## DMRcate, 1045

find_overlaps(bumphunterComp, ProbeLassoComp) # 0, 1_2
find_overlaps(bumphunterComp, DMRcateComp) #0, 1_3
find_overlaps(ProbeLassoComp, DMRcateComp) #67 , 2_3

find_overlaps(find_overlaps(bumphunterComp, ProbeLassoComp), DMRcateComp) #0

Venn_comp <- draw.triple.venn(area1 = 2,
                              area2 = 122,
                              area3 = 1045,
                              n12 = 0,
                              n13 = 0,
                              n23 = 67,
                              n123 = 0,
                              category = c("Bumphunter", "Probe Lasso", "DMRcate"),
                              fill = c("pink", "lightblue", "lightgoldenrodyellow"),
                              cex = 1.2, cat.cex = 0.6, 
                              cat.fontface = 2,
                              cat.just =list(c(0, -1.5), c(0, 1), c(0.7 , 2)),
                              # lty = "blank"
                              cat.dist = c(0.05, 0.05, 0.05)) 

Venn_comp <- grid.arrange(gTree(children=Venn_comp), 
                          top=textGrob("Comparable parameters",
                                        gp=gpar(fontface="bold",cex=0.8))) %>% 
  as_ggplot()

# jpeg(file = here("Method_DNAme/figures/test.jpeg"),
#      width = 8, height = 8, units = "cm", res = 320)
# Venn_comp
# dev.off()



VennPlots_fig <- ggarrange(Venn_def, Venn_opti, Venn_comp,
         labels = c("A", "B", "C", NULL),
         font.label = list(size = 12, face = "plain", color ="black"), 
         nrow=2, ncol = 2)

jpeg(file = here("Method_DNAme/figures/VennPlots_DMRs_placenta.jpeg"),
     width = 16, height = 16, units = "cm", res = 320)
VennPlots_fig
dev.off()


