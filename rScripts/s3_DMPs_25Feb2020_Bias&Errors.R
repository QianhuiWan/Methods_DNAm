# library(doParallel)
# registerDoParallel(cores = 3)

library(here)
library(minfi)
library(tidyverse)
library(magrittr)
library(microbenchmark)
library(missMethyl)
library(limma)
library(gamlss)
# library(betareg)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# BM_list <- readRDS(file = here("Method_DNAme/rds/preprocessed_PAC12_placenta.rds"))

M <- readRDS(file = here("Method_DNAme/rds/Mval.rds"))
B <- readRDS(file = here("Method_DNAme/rds/Bval.rds"))

set.seed(2030)
M <- M[sample(nrow(M), 100),]
B <- B[match(rownames(M), rownames(B)),]

# saveRDS(M, file = here("Method_DNAme/rds/Mval.rds"))
# saveRDS(B, file = here("Method_DNAme/rds/Bval.rds"))

baseDir <- file.path(here("EPICdata/IDAT_Placenta_PAC_sub12"))
phenoData <- read.metharray.sheet(baseDir)
EpicRGsetPAC12 <- read.metharray.exp(targets = phenoData, extended = TRUE, force=TRUE)

# estSex <- EpicRGsetPAC12 %>% preprocessRaw() %>% mapToGenome() %>% getSex()
# phenoDataPAC12 <- EpicRGsetPAC12 %>% preprocessRaw() %>% mapToGenome() %>% 
#   addSex(sex = estSex) %>% 
#   pData()

phenoDataPAC12 <- readRDS(file = here("Method_DNAme/rds/phenoDataPAC12.rds"))

# DMPS with different methods #########################################################

# variables
Trimester <- factor(phenoDataPAC12$Trimester)
GA <- phenoDataPAC12$Gestational.Age %>% as.numeric()
# FetalSex <- factor(phenoDataPAC12$Fetal.Sex, levels = c("F", "M"))
ArrayDateBatch <- factor(phenoDataPAC12$ArrayDateReceived)
## design the matrix for limma
### create a matrix
design <- model.matrix(~ Trimester+ArrayDateBatch)
### rename the col of the matrix
colnames(design) <- c("Intercept","secondTrimester",
                      "arrayDateReceivedBatch")


############## Use functions ############################################################

table(colnames(B)==rownames(phenoDataPAC12))

library(boot)

set.seed(2077)

dmp_lmFUN <- function(dat,indicies){
  dat <- dat[indicies,]
  fit_lm <- eval(bquote(lm.fit(x = design, y = dat))) #replace this line with lm()
  fit_coef <- coef(fit_lm)[2]
  return(fit_coef)
}


dmp_betaFUN <- function(dat,indicies){
  dat <- dat[indicies]
  fit_beta <- eval(bquote(gamlss(formula = dat~phenoDataPAC12$Trimester+phenoDataPAC12$ArrayDateReceived, 
                                 family = BE(), 
                                 control = gamlss.control(trace = FALSE)))) #replace this line with betareg (or gamlss)
  fit_coef <- coef(fit_beta)[2]
  return(fit_coef)
}


lm_boot <- boot(data=t(M),statistic=dmp_lmFUN,R=500,parallel="multicore",ncpus=3)

glmBeta_boot <- boot(data=t(B),statistic=dmp_betaFUN,R=500,parallel="multicore",ncpus=3)

# The bootstrap summary:
lm_boot$t
glmBeta_boot$t
# To get the confidence intervals:
boot.ci(lm_boot,type="basic") #LM coef
boot.ci(glmBeta_boot,type="basic") #RUV coef


