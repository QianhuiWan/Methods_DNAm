## use multiple cores
library(doParallel)
registerDoParallel(cores = 3)

library(here)
library(rngtools)
library(minfi)
library(tidyverse)
library(magrittr)
library(microbenchmark)
library(limma)
library(ChAMP)
library(DMRcate)
library(plyranges)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(pryr)

BM_list <- readRDS(file = here("Method_DNAme/rds/preprocessed_GEO4.rds"))

M <- BM_list$Mnorm
B <- BM_list$betaNorm

phenoData <- BM_list$phenoDataGEO

# DMRs with different methods #########################################################
# variables
Cell <- str_replace_all(phenoData$`Cell line`, c("WI-38 fibroblasts expressing control"="control", "WI-38 fibroblasts overexpressing NNMT"="NNMToverEx")) %>% 
  factor(., levels = c("control", "NNMToverEx"))

## design the matrix for limma
### create a matrix
design <- model.matrix(~ Cell)
### rename the col of the matrix
colnames(design) <- c("Intercept","NNMToverEx")

############ bumphunter ############################################################
# make genomic ratio set for minfi::bumphunter
bumphunter_DMR <- function(B, design){

GRset <- makeGenomicRatioSetFromMatrix(B, array = "IlluminaHumanMethylationEPIC",
                                          annotation = "ilm10b4.hg19", mergeManifest = TRUE,
                                          what = "Beta")

## hunt our DMRs, linear, so we should use beta values (bumphunter will convert it into M values)
### minfi bumphunter
bump_dmrs <- minfi::bumphunter(GRset, design=design,
                               coef=2,
                               # cutoff=0.2,
                               pickCutoff=TRUE,
                               pickCutoffQ = 0.99, # this is for permutation
                               nullMethod="bootstrap",
                               smooth=TRUE,
                               B=5,
                               maxGap=500,
                               smoothFunction=loessByCluster,
                               minNum=2,
                               useWeights=TRUE)

return(bump_dmrs)
}

######### Probe Lasso ############################################################################

ProbeLasso_DMR <- function(B, phenoData){
  
Cell <- str_replace_all(phenoData$`Cell line`, c("WI-38 fibroblasts expressing control"="control", "WI-38 fibroblasts overexpressing NNMT"="NNMToverEx")) %>% 
       factor(., levels = c("control", "NNMToverEx"))

### ProbeLasso found 108 DMRs with P<=0.05
probeLasso_dmrs <- champ.DMR(beta = B,
          pheno = Cell,
          adjPvalProbe = 0.05,
          arraytype = "EPIC",
          method = "ProbeLasso",
          minDmrSize = 50,
          minDmrSep =500,
          minProbes = 2,
          meanLassoRadius=1000)

return(probeLasso_dmrs)
}

######### DMRcate ################################################################################

DMRcate_DMR <- function(M, design){
  
colnames(design) <- c("(Intercept)","NNMToverEx")

# estimate weights for each sample
w <- arrayWeights(M, design = design)

# Create genomic ratio sets
GRset <- makeGenomicRatioSetFromMatrix(M, array = "IlluminaHumanMethylationEPIC",
                                       annotation = "ilm10b4.hg19", mergeManifest = TRUE,
                                       what = "M")

# cpg.annotateFUN
cpgAnnotation <- cpg.annotate(object = GRset,
                              datatype = "array",
                              what = "M",
                              analysis.type = "differential",
                              design = design,
                              # contrasts = TRUE,
                              # cont.matrix = contMatrix,
                              fdr = 0.05,
                              adjust.method ="BH",
                              coef = "NNMToverEx",
                              weights = w,
                              arraytype = "EPIC")

## run dmrcate with multi-core(cores=5)
DMRcate_dmrs <- dmrcate(object = cpgAnnotation,
                              lambda = 500,
                              C = 5,
                              min.cpgs = 2)

return(DMRcate_dmrs)
}

############## Use functions ############################################################

t1 <- system.time(bumphunter <-  bumphunter_DMR(B = B, design = design))
t2 <- system.time(ProbeLasso <-  ProbeLasso_DMR(B = B, phenoData = phenoData))
t3 <- system.time(DMRcate <-  DMRcate_DMR(M = M, design = design))

mem_used()
object_size(bumphunter)
object_size(ProbeLasso)
object_size(DMRcate)

save(bumphunter, ProbeLasso, DMRcate,
     file = here("Method_DNAme/rds/DMRmethods_ComparablePara_GEO.RData"))
# 
# mbm <- microbenchmark(
#    bumphunter =  bumphunter_DMR(B = B, design = design),
#    ProbeLasso =  ProbeLasso_DMR(B = B, phenoData = phenoData),
#    DMRcate =  DMRcate_DMR(M = M, design = design), times = 5
#  )
# 
# library(ggplot2)
# library(ggpubr)
# theme_set(theme_pubr())
# mbm_plot <- autoplot(mbm)

# ggsave(filename = "mbmPlot_DMRs.jpeg",
#         plot = mbm_plot,
#         device = "jpeg",
#         path = here("Method_DNAme/rds_figure/"),scale = 1.6 ,
#         width = 18, height = 5.5, units = "cm", dpi = 320)

# saveRDS(mbm_plot, file=here("Method_DNAme/rds_figure/mbm_GEO_DMRcomparablePara_plot.rds"))


# mbm_plot_GEOdefault <- readRDS(file=here("Method_DNAme/rds_figure/mbm_GEO_DMR_defaultPara_plot.rds"))+
#     scale_y_continuous(limits = c(0, 3500), name = "Time [seconds]")+coord_flip()+geom_boxplot()
# 
# mbm_plot_GEOopti <- readRDS(file=here("Method_DNAme/rds_figure/mbm_GEO_DMRoptiPara_plot.rds"))+
#     scale_y_continuous(limits = c(0, 3500), name = "Time [seconds]")+coord_flip()+geom_boxplot()
# 
# mbm_plot_GEOcomp <- readRDS(file=here("Method_DNAme/rds_figure/mbm_GEO_DMRcomparablePara_plot.rds"))+
#     scale_y_continuous(limits = c(0, 3500), name = "Time [seconds]")+coord_flip()+geom_boxplot()
# 
# 
# Figure4_method <- ggarrange(mbm_plot_GEOdefault, mbm_plot_GEOopti, mbm_plot_GEOcomp,
#                             labels = c("A", "B", "C"),
#                             font.label = list(size = 12, face = "plain", color ="black"),
#                             nrow = 3, ncol=1)
# 
# jpeg(file = here("Method_DNAme/figures/Figure4_method_v1.jpeg"),
#      width = 17, height = 12, units = "cm", res = 320)
# Figure4_method
# 
# dev.off()

########### true DMRs #######################################
# bumphunter, ProbeLasso, DMRcate
load(file = here("Method_DNAme/rds/DMRmethods_ComparablePara_GEO.RData"))

## bumphunter
bumphunterComp <- bumphunter$table[bumphunter$table$fwer<0.05,] %>% 
  dplyr::select(seqnames=chr, start, end, TotalProbeNum=L, deltaB=value) %>% 
  magrittr::inset("index", value=paste0("DMR", 1:nrow(.))) %>% 
  as_granges() 
table(abs(bumphunterComp$deltaB) >0.2)

## ProbeLasso
ProbeLassoComp <- ProbeLasso$ProbeLassoDMR %>% tibble::rownames_to_column(var = "index") %>% 
  as_granges() %>% 
  plyranges::filter(dmrP<0.05)

table(abs(ProbeLassoComp$betaAv_control-ProbeLassoComp$betaAv_NNMToverEx) >0.2)

## DMRcate
DMRcateComp <- extractRanges(DMRcate, genome = "hg19") %>% 
  plyranges::filter(Stouffer<0.05)
DMRcateComp$index <- paste0("DMR", 1:length(DMRcateComp))

table(abs(DMRcateComp$meandiff) >0.2)



########## overlapping with promoters and enhancers ########################
## GR of promoters (184476) and enhancers (32693)
Promoter_humanAll <- read.csv(here("PrepareData/DownloadedFilesForEPICpipeline", "promoter_data_at_2018-07-13_10-56-38_allHuman.bed"), sep = "\t", stringsAsFactors = FALSE,header = FALSE,
                              col.names = c("seqnames", "start", "end","coord", "V5", "V6")) %>% as_granges()

Enhancer_humanAll <- read.csv(here("PrepareData/DownloadedFilesForEPICpipeline", "enhancer_data_at_2018-07-13_10-17-59_allHuman.bed"), sep = "\t",stringsAsFactors = FALSE,header = FALSE,
                              col.names = c("seqnames", "start", "end","coord", "V5", "V6")) %>% as_granges()

## bumphunter
bumphunterComp <- bumphunterComp %>% plyranges::filter(abs(deltaB)>0.2)
find_overlaps_within(Promoter_humanAll, bumphunterComp)$index %>% unique() %>% length()
find_overlaps_within(Enhancer_humanAll, bumphunterComp)$index %>% unique() %>% length()

## ProbeLasso
ProbeLassoComp <- ProbeLassoComp %>% 
  plyranges::filter(abs(ProbeLassoComp$betaAv_control-ProbeLassoComp$betaAv_NNMToverEx) >0.2)
find_overlaps_within(Promoter_humanAll, ProbeLassoComp)$index %>% unique() %>% length()
find_overlaps_within(Enhancer_humanAll, ProbeLassoComp)$index %>% unique() %>% length()

# DMRcate
DMRcateComp <- DMRcateComp %>% plyranges::filter(abs(meandiff) >0.2)
find_overlaps_within(Promoter_humanAll, DMRcateComp)$index %>% unique() %>% length()
find_overlaps_within(Enhancer_humanAll, DMRcateComp)$index %>% unique() %>% length()

save(bumphunterComp, ProbeLassoComp, DMRcateComp, 
     file = here("Method_DNAme/rds/DMRmethods_ComparablePara_GEO_Grange.RData"))


