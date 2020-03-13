
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
# bumphunter, ProbeLasso, DMRcate
load(file = here("Method_DNAme/rds/DMRs_methods.RData"))

# bumphunter_def, ProbeLasso_def, DMRcate_def,
load(file = here("Method_DNAme/rds/DMRs_defaultPara_methods.RData"))

# optimised para
## bumphunter, 64
bumphunter_optiGR <- bumphunter$table %>% 
  dplyr::select(seqnames=chr, start, end, TotalProbeNum=L, deltaB=value) %>% 
  as_granges() 

## Probe Lasso, 2
ProbeLasso_optiGR <- ProbeLasso$ProbeLassoDMR %>% 
  as_granges() %>% 
  plyranges::select(dmrP, dmrSize) 

## DMRcate, 260
DMRcate_optiGR <- extractRanges(DMRcate, genome = "hg19") %>% 
  plyranges::filter(Stouffer<0.05&abs(meandiff)>0.2) %>% 
  plyranges::select(no.cpgs, Stouffer, meandiff) 

find_overlaps(bumphunter_optiGR, ProbeLasso_optiGR) #1, 1_2
find_overlaps(bumphunter_optiGR, DMRcate_optiGR) #12, 1_3
find_overlaps(ProbeLasso_optiGR, DMRcate_optiGR) #0, 2_3

find_overlaps(find_overlaps(bumphunter_optiGR, ProbeLasso_optiGR), DMRcate_optiGR) #0

Venn_opti <- draw.triple.venn(area1 = 64,
                              area2 = 2,
                              area3 = 260,
                              n12 = 1,
                              n13 = 12,
                              n23 = 0,
                              n123 = 0,
                              category = c("Bumphunter", "Probe Lasso", "DMRcate"),
                              fill = c("pink", "lightblue", "lightgoldenrodyellow"),
                              cex = 1.5, cat.cex = 1, 
                              cat.fontface = 2,
                              cat.just =list(c(1.2, -1), c(0.4, -0.8), c(0, -1)),
                              # lty = "blank"
                              cat.dist = c(0.05, -0.4, 0.05)) 

Venn_opti <- grid.arrange(gTree(children=Venn_opti), 
                          top=textGrob("Optimised parameters",
                                        gp=gpar(fontface="bold",cex=1))) %>% 
  as_ggplot()



##############################################################################
# default para
## bumphunter, 209
bumphunter_def_GR <- bumphunter_def$table %>% 
  dplyr::select(seqnames=chr, start, end, TotalProbeNum=L, deltaB=value) %>% 
  as_granges() 

## Probe Lasso, 1
ProbeLasso_def_GR <- ProbeLasso_def$ProbeLassoDMR %>% 
  as_granges() %>% 
  plyranges::select(dmrP, dmrSize) 


## DMRcate, 369
DMRcate_def_GR <- extractRanges(DMRcate_def, genome = "hg19") %>% 
  plyranges::filter(Stouffer<0.05&abs(meandiff)>0.2) %>% 
  plyranges::select(no.cpgs, Stouffer, meandiff) 

find_overlaps(bumphunter_def_GR, ProbeLasso_def_GR) #0, 1_2
find_overlaps(bumphunter_def_GR, DMRcate_def_GR) #14, 1_3
find_overlaps(ProbeLasso_def_GR, DMRcate_def_GR) #0, 2_3

find_overlaps(find_overlaps(bumphunter_def_GR, ProbeLasso_def_GR), DMRcate_def_GR) #0


# Venn diagram
Venn_def <- draw.triple.venn(area1 = 209,
                              area2 = 1,
                              area3 = 369,
                              n12 = 0,
                              n13 = 14,
                              n23 = 0,
                              n123 = 0,
                              category = c("Bumphunter", "Probe Lasso", "DMRcate"),
                              fill = c("pink", "lightblue", "lightgoldenrodyellow"),
                              cex = 1.3, cat.cex = 1, 
                              cat.fontface = 2,
                              cat.just =list(c(1.2, -1), c(0.5, 1), c(0, -1)),
                              # lty = "blank"
                              cat.dist = c(0.05, 0.05, 0.05)
                             ) 

Venn_def <- grid.arrange(gTree(children=Venn_def), 
                          top=textGrob("Default parameters",
                                       gp=gpar(fontface="bold",cex=1))) %>% 
  as_ggplot()


VennPlots_fig <- ggarrange(Venn_def, Venn_opti,
         labels = c("A", "B"),
         font.label = list(size = 12, face = "plain", color ="black"), 
         nrow=1, ncol = 2)

jpeg(file = here("Method_DNAme/figures/VennPlots_DMRs.jpeg"),
     width = 16, height = 8, units = "cm", res = 320)
VennPlots_fig
dev.off()


