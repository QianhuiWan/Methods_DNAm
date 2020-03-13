## use multiple cores in freya
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
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(pryr)

# BM_list <- readRDS(file = here("Method_DNAme/rds/preprocessed_PAC12_placenta.rds"))

M <- readRDS(file = here("Method_DNAme/rds/Mval.rds"))
B <- readRDS(file = here("Method_DNAme/rds/Bval.rds"))


phenoDataPAC12 <- readRDS(file = here("Method_DNAme/rds/phenoDataPAC12.rds"))

# DMRs with different methods #########################################################
# variables
Trimester <- factor(phenoDataPAC12$Trimester)
# FetalSex <- factor(phenoDataPAC12$Fetal.Sex, levels = c("F", "M"))
ArrayDateBatch <- factor(phenoDataPAC12$ArrayDateReceived)
## design the matrix for limma
### create a matrix
design <- model.matrix(~ Trimester+ArrayDateBatch)
### rename the col of the matrix
colnames(design) <- c("Intercept","secondTrimester",
                      "arrayDateReceivedBatch")

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
                               cutoff=0.2,
                               nullMethod="bootstrap",
                               smooth=TRUE,
                               # pickCutoffQ = 0.95, # this is for permutation
                               maxGap=500,
                               smoothFunction=loessByCluster,
                               useWeights=TRUE)

return(bump_dmrs)
}

######### Probe Lasso ############################################################################

ProbeLasso_DMR <- function(B, phenoDataPAC12){

### ProbeLasso found 108 DMRs with P<=0.05
probeLasso_dmrs <- champ.DMR(beta = B,
          pheno = phenoDataPAC12$Trimester,
          # adjPvalProbe = 0.05,
          arraytype = "EPIC",
          method = "ProbeLasso",
          # minDmrSize = 50,
          # minDmrSep =1000,
          # meanLassoRadius=1000,
          cores = 3)

return(probeLasso_dmrs)
}

######### DMRcate ################################################################################

DMRcate_DMR <- function(M, design){
  
colnames(design) <- c("(Intercept)","secondTrimester",
                      "arrayDateReceivedBatch")

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
                              coef = "secondTrimester",
                              weights = w,
                              arraytype = "EPIC")

## run dmrcate with multi-core(cores=5)
DMRcate_dmrs <- dmrcate(object = cpgAnnotation,
                              # lambda = 500,
                              # C = 5,
                              # min.cpgs = 5,
                              min.cpgs = 3)

return(DMRcate_dmrs)
}

############## Use functions ############################################################

t1 <- system.time(bumphunter_def <-  bumphunter_DMR(B = B, design = design))
t2 <- system.time(ProbeLasso_def <-  ProbeLasso_DMR(B = B, phenoDataPAC12 = phenoDataPAC12))
t3 <- system.time(DMRcate_def <-  DMRcate_DMR(M = M, design = design))



mem_used()
object_size(bumphunter_def)
object_size(ProbeLasso_def)
object_size(DMRcate_def)

save(bumphunter_def, ProbeLasso_def, DMRcate_def, 
     file = here("Method_DNAme/rds/DMRs_defaultPara_methods.RData"))

# mbm_default <- microbenchmark(
#   bumphunter =  bumphunter_DMR(B = B, design = design),
#   ProbeLasso =  ProbeLasso_DMR(B = B, phenoDataPAC12 = phenoDataPAC12),
#   DMRcate =  DMRcate_DMR(M = M, design = design), times = 20
# )
# 
# library(ggplot2)
# library(ggpubr)
# theme_set(theme_pubr())
# mbm_default_plot <- autoplot(mbm_default)
# 
# ggsave(filename = "mbmPlot_DMRs_default.jpeg",
#        plot = mbm_default_plot,
#        device = "jpeg",
#        path = here("Method_DNAme/rds_figure/"),scale = 1.6 ,
#        width = 18, height = 5.5, units = "cm", dpi = 320)

# saveRDS(mbm_default_plot, file = here("Method_DNAme/rds_figure/mbm_DMRdefault_plot.rds"))


