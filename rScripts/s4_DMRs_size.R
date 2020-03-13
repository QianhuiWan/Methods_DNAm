
library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(plyranges)
library(minfi)
# bumphunter, ProbeLasso, DMRcate
load(file = here("Method_DNAme/rds/DMRs_methods.RData"))

# bumphunter_def, ProbeLasso_def, DMRcate_def,
load(file = here("Method_DNAme/rds/DMRs_defaultPara_methods.RData"))

# Get Probe GR object
B <- readRDS(file = here("Method_DNAme/rds/Bval.rds"))

GRset <- makeGenomicRatioSetFromMatrix(B, array = "IlluminaHumanMethylationEPIC",
                                       annotation = "ilm10b4.hg19", mergeManifest = TRUE,
                                       what = "Beta")

B_ProbeGR <- GRset@rowRanges[, c(1:2,7)]

matchWithGRprobe <- base::match(B_ProbeGR$Name, rownames(B))
values(B_ProbeGR) <- cbind(values(B_ProbeGR), DataFrame(B[matchWithGRprobe,], check.names = F)) 


# optimised para
## bumphunter, 64
bumphunter_opti <- bumphunter$table %>% 
  dplyr::select(seqnames=chr, start, end, TotalProbeNum=L, deltaB=value) %>% 
  as_granges() %>% 
  group_by_overlaps(B_ProbeGR[, "Name"]) %>% 
  plyranges::summarise(RealProbeNum=n()) %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("Bumphunter", nrow(.)))

## Probe Lasso, 2
ProbeLasso_opti <- ProbeLasso$ProbeLassoDMR %>% 
  as_granges() %>% 
  group_by_overlaps(B_ProbeGR[, "Name"]) %>% 
  plyranges::summarise(RealProbeNum=n()) %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("Probe Lasso", nrow(.)))

## DMRcate, 2111
DMRcate_opti <- extractRanges(DMRcate, genome = "hg19") %>% 
  plyranges::filter(Stouffer<0.05&abs(meandiff)>0.2) %>% 
  plyranges::select(no.cpgs, Stouffer, meandiff) %>% 
  group_by_overlaps(B_ProbeGR[, "Name"]) %>% 
  plyranges::summarise(RealProbeNum=n()) %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("DMRcate", nrow(.)))

Three_opti <- rbind(bumphunter_opti, ProbeLasso_opti, DMRcate_opti) %>% 
  magrittr::inset("Parameter", value=rep("Optimised", nrow(.)))

# default para
## bumphunter, 209
bumphunter_def_PN <- bumphunter_def$table %>% 
  dplyr::select(seqnames=chr, start, end, TotalProbeNum=L, deltaB=value) %>% 
  as_granges() %>% 
  group_by_overlaps(B_ProbeGR[, "Name"]) %>% 
  plyranges::summarise(RealProbeNum=n()) %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("Bumphunter", nrow(.)))

## Probe Lasso, 1
ProbeLasso_def_PN <- ProbeLasso_def$ProbeLassoDMR %>% 
  as_granges() %>% 
  group_by_overlaps(B_ProbeGR[, "Name"]) %>% 
  plyranges::summarise(RealProbeNum=n()) %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("Probe Lasso", nrow(.)))


## DMRcate, 3410
DMRcate_def_PN <- extractRanges(DMRcate_def, genome = "hg19") %>% 
  plyranges::filter(Stouffer<0.05&abs(meandiff)>0.2) %>% 
  plyranges::select(no.cpgs, Stouffer, meandiff) %>% 
  group_by_overlaps(B_ProbeGR[, "Name"]) %>% 
  plyranges::summarise(RealProbeNum=n()) %>% as.data.frame() %>% 
  magrittr::inset("Methods", value=rep("DMRcate", nrow(.)))

Three_def <- rbind(bumphunter_def_PN, ProbeLasso_def_PN, DMRcate_def_PN) %>% 
  magrittr::inset("Parameter", value=rep("Default", nrow(.)))

# rbind all rows together
theme_set(theme_pubr())
DMR_size_plot <- rbind(Three_opti,Three_def) %>% 
  ggplot(aes(x=Methods, y=RealProbeNum, fill=Parameter))+
  geom_boxplot()+
  # facet_grid(.~Methods)+
  # geom_violin()+
  stat_compare_means(label = "p.signif", method = "wilcox.test", size = 5)+
  # geom_signif(comparisons = list(c("Default", "Optimised")), map_signif_level=TRUE)+ #,textsize = 6 , test = "wilcox.test"
  labs(y="Number of cites in DMRs")+
  theme(
    text = element_text(size=12, colour = "black"),
    axis.text.x = element_text(size=12, colour = "black"),
    axis.title.x = element_text(size=12, colour = "black"),
    axis.text.y=element_text(size=12, colour = "black"), 
    legend.position = "right")

ggsave(filename = "DMR_size_default&Opti_v2.jpeg",
       plot = DMR_size_plot,
       device = "jpeg",
       path = here("Method_DNAme/rds_figure/"),scale = 1.6 ,
       width = 17, height = 6, units = "cm", dpi = 320)


