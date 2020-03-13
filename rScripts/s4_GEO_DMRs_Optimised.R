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
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(pryr)

BM_list <- readRDS(file = here("Method_DNAme/rds/preprocessed_GEO4.rds"))

M <- BM_list$Mnorm
B <- BM_list$betaNorm

phenoData <- BM_list$phenoDataGEO

# DMRs with different methods #########################################################
# variables
Cell <- factor(phenoData$`Cell line`)

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
                               pickCutoffQ = 0.95, # this is for permutation
                               nullMethod="bootstrap",
                               smooth=TRUE,
                               B=5,
                               maxGap=250,
                               smoothFunction=loessByCluster,
                               useWeights=TRUE)

return(bump_dmrs)
}

######### Probe Lasso ############################################################################

ProbeLasso_DMR <- function(B, phenoData){

### ProbeLasso found 108 DMRs with P<=0.05
probeLasso_dmrs <- champ.DMR(beta = B,
          pheno = phenoData$`Cell line`,
          adjPvalProbe = 0.05,
          arraytype = "EPIC",
          method = "ProbeLasso",
          minDmrSize = 50,
          minDmrSep =1000,
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
     file = here("Method_DNAme/rds/DMRmethods_optiPara_GEO.RData"))

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
#  
# ggsave(filename = "mbmPlot_DMRs.jpeg",
#         plot = mbm_plot,
#         device = "jpeg",
#         path = here("Method_DNAme/rds_figure/"),scale = 1.6 ,
#         width = 18, height = 5.5, units = "cm", dpi = 320)
# 
# saveRDS(mbm_plot, file=here("Method_DNAme/rds_figure/mbm_DMRoptiPara_plot.rds"))

# mbm_plot <- readRDS(file=here("Method_DNAme/rds_figure/mbm_DMRoptiPara_plot.rds"))+
#   scale_y_continuous(limits = c(0, 200), name = "Time [seconds]")+coord_flip()
# 
# mbm_default_plot <- readRDS(file = here("Method_DNAme/rds_figure/mbm_DMRdefault_plot.rds"))+
#   scale_y_continuous(limits = c(0, 200), name = "Time [seconds]")+coord_flip()
# 
# 
# Figure4_method <- ggarrange(mbm_default_plot, mbm_plot, labels = c("A", "B"), 
#                             font.label = list(size = 12, face = "plain", color ="black"), 
#                             nrow = 2, ncol=1)
# 
# jpeg(file = here("Method_DNAme/figures/Figure4_method.jpeg"),
#      width = 17, height = 12, units = "cm", res = 320)
# Figure4_method
# 
# dev.off()
