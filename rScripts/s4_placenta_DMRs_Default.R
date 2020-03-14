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

BM_list <- readRDS(file = here("Method_DNAme/rds/preprocessed_PAC12_placenta.rds"))

M <- BM_list$Mnorm
B <- BM_list$betaNorm
phenoData <- BM_list$phenoDataPAC12

# DMRs with different methods #########################################################
# variables
Trimester <- factor(phenoData$Trimester)
ArrayDateBatch <- factor(phenoData$ArrayDateReceived)

## design the matrix for limma
### create a matrix
design <- model.matrix(~ Trimester+ArrayDateBatch)
### rename the col of the matrix
colnames(design) <- c("Intercept","T2", "Batch")

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
                               # smooth=TRUE,
                               B=5,
                               maxGap=500,
                               smoothFunction=loessByCluster,
                               minNum=7,
                               useWeights=FALSE)

return(bump_dmrs)
}

######### Probe Lasso ############################################################################

ProbeLasso_DMR <- function(B, phenoData){

### ProbeLasso found 108 DMRs with P<=0.05
probeLasso_dmrs <- champ.DMR(beta = B,
          pheno = phenoData$Trimester,
          adjPvalProbe = 0.05,
          arraytype = "EPIC",
          method = "ProbeLasso",
          minDmrSize = 50,
          minDmrSep =1000,
          minProbes = 7,
          meanLassoRadius=375)

return(probeLasso_dmrs)
}

######### DMRcate ################################################################################

DMRcate_DMR <- function(M, design){

colnames(design) <- c("(Intercept)","T2", "Batch")

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
                              coef = "T2",
                              weights = w,
                              arraytype = "EPIC")

## run dmrcate with multi-core(cores=5)
DMRcate_dmrs <- dmrcate(object = cpgAnnotation,
                              lambda = 1000,
                              C = 2,
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
     file = here("Method_DNAme/rds/DMRmethods_DefaultPara.RData"))

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


# saveRDS(mbm_plot, file=here("Method_DNAme/rds_figure/mbm_GEO_DMR_defaultPara_plot.rds"))

########### true DMRs #######################################
# bumphunter, ProbeLasso, DMRcate, 
load(file = here("Method_DNAme/rds/DMRmethods_DefaultPara.RData"))

## bumphunter
bumphunterDefault <- bumphunter$table[bumphunter$table$fwer<0.05,] %>% 
  dplyr::select(seqnames=chr, start, end, TotalProbeNum=L, deltaB=value) %>% 
  magrittr::inset("index", value=paste0("DMR", 1:nrow(.))) %>% 
  as_granges() 
table(abs(bumphunterDefault$deltaB) >0.2)

## ProbeLasso
ProbeLassoDefault <- ProbeLasso$ProbeLassoDMR %>% tibble::rownames_to_column(var = "index") %>% 
  as_granges() %>% 
  plyranges::filter(dmrP<0.05)
  
table(abs(ProbeLassoDefault$betaAv_Second-ProbeLassoDefault$betaAv_First) >0.2)

## DMRcate
DMRcateDefault <- extractRanges(DMRcate, genome = "hg19") %>% 
  plyranges::filter(Stouffer<0.05)
DMRcateDefault$index <- paste0("DMR", 1:length(DMRcateDefault))

table(abs(DMRcateDefault$meandiff) >0.2)



########## overlapping with promoters and enhancers ########################
## GR of promoters (184476) and enhancers (32693)
Promoter_humanAll <- read.csv(here("PrepareData/DownloadedFilesForEPICpipeline", "promoter_data_at_2018-07-13_10-56-38_allHuman.bed"), sep = "\t", stringsAsFactors = FALSE,header = FALSE,
                              col.names = c("seqnames", "start", "end","coord", "V5", "V6")) %>% as_granges()

Enhancer_humanAll <- read.csv(here("PrepareData/DownloadedFilesForEPICpipeline", "enhancer_data_at_2018-07-13_10-17-59_allHuman.bed"), sep = "\t",stringsAsFactors = FALSE,header = FALSE,
                              col.names = c("seqnames", "start", "end","coord", "V5", "V6")) %>% as_granges()

## bumphunter
bumphunterDefault <- bumphunterDefault %>% plyranges::filter(abs(deltaB)>0.2)
find_overlaps_within(Promoter_humanAll, bumphunterDefault)$index %>% unique() %>% length()
find_overlaps_within(Enhancer_humanAll, bumphunterDefault)$index %>% unique() %>% length()

## ProbeLasso
ProbeLassoDefault <- ProbeLassoDefault %>% 
  plyranges::filter(abs(ProbeLassoDefault$betaAv_Second-ProbeLassoDefault$betaAv_First) >0.2)
find_overlaps_within(Promoter_humanAll, ProbeLassoDefault)$index %>% unique() %>% length()
find_overlaps_within(Enhancer_humanAll, ProbeLassoDefault)$index %>% unique() %>% length()

# DMRcate
DMRcateDefault <- DMRcateDefault %>% plyranges::filter(abs(meandiff) >0.2)
find_overlaps_within(Promoter_humanAll, DMRcateDefault)$index %>% unique() %>% length()
find_overlaps_within(Enhancer_humanAll, DMRcateDefault)$index %>% unique() %>% length()


save(bumphunterDefault, ProbeLassoDefault, DMRcateDefault, 
     file = here("Method_DNAme/rds/DMRmethods_DefaultPara_Grange.RData"))

