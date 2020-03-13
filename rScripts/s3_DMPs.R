
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

BM_list <- readRDS(file = here("Method_DNAme/rds/preprocessed_PAC12_placenta.rds"))

baseDir <- file.path(here("EPICdata/IDAT_Placenta_PAC_sub12"))
phenoData <- read.metharray.sheet(baseDir)
EpicRGsetPAC12 <- read.metharray.exp(targets = phenoData, extended = TRUE, force=TRUE)

estSex <- EpicRGsetPAC12 %>% preprocessRaw() %>% mapToGenome() %>% getSex()
phenoDataPAC12 <- EpicRGsetPAC12 %>% preprocessRaw() %>% mapToGenome() %>% 
  addSex(sex = estSex) %>% 
  pData()

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

############ RUVm ############################################################
# RUVm Illumina negtive controls (INCs): ====
## true negtive control probes, 411
## INCs are Illumina negtive controls # M=BM_list$Mnorm

RUVm_dmp <- function(EpicRGsetPAC12, M){

INCs <- getINCs(EpicRGsetPAC12)

## Mc means combine m-value and Illumina negtive controls
Mc <- base::rbind(M, INCs)

## ctl means ctl probes' information
ctl1 <- base::rownames(Mc) %in% base::rownames(INCs)

# RUVfit and RUVadj====
#fit microarray linear regression model
rfit1 <- missMethyl::RUVfit(Y = Mc, X= Trimester, ctl = ctl1)

# empirical bayes method: adjust variance
rfit2 <- RUVadj(Y = Mc, fit = rfit1)

#Extract a table of the top-ranked CpGs from a linear model 
top1 <- topRUV(rfit2, number = Inf, p.BH = 1)
# head(top1)

## ECPs, Empirical control probes
### p.ebayes.BH is FDR-adjusted p-values, after applying the empirical bayes method

# probes without significant differenc as EPCs====
ctl_EPCs <- base::rownames(M) %in% base::rownames(top1[top1$p.BH_X1.Second >0.5,])
# table(ctl_EPCs)

# use ECPs to do 2nd Diff methylation====
# fit linear model, 2nd time
rfit1_2nd <- RUVfit(Y = M, X= Trimester, ctl = ctl_EPCs)

# adjust variance
rfit2_2nd <- RUVadj(Y = M, rfit1_2nd)

# extract topRUV result====
top2 <- topRUV(rfit2_2nd, number = Inf)  

table(top2$p.BH_X1.Second <0.05) # 35937

return(top2)

}

######### limma ################################################################
Lm_dmp <- function(phenoDataPAC12, M){
# design the matrix for limma

# estimate weights for each sample
w <- arrayWeights(M, design = design)

phenoDataPAC12 %>%
  cbind(w) %>%
  as.data.frame() %>%
  as_tibble() %>%
  ggplot(aes(Trimester, w)) +
  geom_violin() +
  geom_jitter(width = 0.2) +
  geom_hline(yintercept = 1, colour = "blue", linetype = 2) +
  facet_wrap(~ArrayDateReceived) +
  theme_bw()

## fit lm using limma: 
fit1_limma <- lmFit(M, design = design, weights = w) %>%
  eBayes()

# decide test
fit1_limma %>%
  decideTests(p.value = 0.05, adjust.method = "BH") %>%
  summary()

## add annotation to top_limma:
annEPIC <- ann850k %>% as.data.frame()
annEPICsub <- annEPIC[match(rownames(M), annEPIC$Name), 
                      c(1:4, 12:19, 22:ncol(annEPIC))]

table(rownames(M)==annEPICsub$Name)

## save BY ajusted data
DMPs_limma <- topTable(fit1_limma, coef = 2,
                       number = Inf, adjust.method="BH",
                       genelist = annEPICsub)

rm(fit1_limma)

return(DMPs_limma)

}

######## Beta regression ########################################################

# Beta regression with `gamlss` R package

# if the meth and unmeth status are individual for each probe, 
# the beta-values obey beta distribution, since beta ∈ {0,1} and 
# beta contain maginal values e.g. round 0 or 1

# use muti cores and also add ajusted covariants in the model:
# library(doMC)
# registerDoMC(cores=3)

##functions

w <- arrayWeights(M, design = design)

## beta glm function
# beta.glm <- function(dat.run, dat.Y, dat.gp, dat.co1){ #, dat.co2
#   dat.tmp <- data.frame("Y"=dat.Y,"X"=dat.gp, "W1"=dat.co1) #, "W2"=dat.co2
#   zero.tab <- dat.tmp$Y<=0.02 #0.02 cutoff
#   zero.prop <- length(zero.tab[zero.tab%in%TRUE])/nrow(dat.tmp)
#   
#   #beta regression
#   #beta
#   if(zero.prop<0.8){fit1 <- gamlss(Y~X+W1,data=dat.tmp,family=BE())} #+W2
#   # if I use betareg package for this fitting, it can be quicker and more comparable with other 2 methods.
#   
#   #zero-inflated
#   if(zero.prop>=0.8){fit1 <- gamlss(Y~X+W1,data=dat.tmp,family=BEZI())}
#   
#   return(fit1)
# }

beta.glm <- function(dat.Y, dat.gp, dat.co1){ #, dat.co2
  # browser()
  dat.tmp <- data.frame("Y"=dat.Y,"X"=dat.gp, "W1"=dat.co1) #, "W2"=dat.co2
  zero.tab <- dat.tmp$Y<=0.02 #0.02 cutoff
  zero.prop <- mean(zero.tab)
  zi <- zero.prop >= 0.8
  fm <- list(BE(), BEZI())[[zi + 1]]
  gamlss(Y~X+W1, data = dat.tmp, family = fm, control = gamlss.control(trace = FALSE))
  # #beta regression
  # #beta
  # if(zero.prop<0.8){fit1 <- gamlss(Y~X+W1,data=dat.tmp,family=BE())} #+W2
  # #zero-inflated
  # if(zero.prop>=0.8){fit1 <- gamlss(Y~X+W1,data=dat.tmp,family=BEZI())}
  # 
  # return(fit1)
}

## calculte zero proportion
CalLowProp <- function(x){
  zero.tab <- x<=0.02
  zero.prop <- length(zero.tab[zero.tab%in%TRUE])/length(x)
  return(zero.prop)
}


## run from here
## we can't use doMC on phoenix, so I need to rewrite 
## this in to a for loop or any other ways that can deal with it

glmBeta_dmp <- function(phenoDataPAC12, B){

dat.run <- B
dat.gp <- model.matrix(~Trimester, data=phenoDataPAC12)[,2] %>% as.matrix()
dat.co1 <- model.matrix(~ArrayDateReceived, data=phenoDataPAC12)[,2] %>% as.matrix()
p.dat <- list()
# dat.co2 <- factor(phenoDataPlacenta$Fetal.Sex, levels = c("F", "M"))

zero.prop <- apply(B, 1, FUN = CalLowProp)

fit1 <- betareg.fit(x = dat.gp, y = dat.run[-which(zero.prop>=0.8), ])

fit1 <- gamlss(B[1,]~Trimester+ArrayDateBatch,family=BE())

for(k in 1:3){ 
  fit_sum <- beta.glm(dat.run[k,],dat.Y = , dat.gp=Trimester, dat.co1=ArrayDateBatch) %>% summary()#, dat.co2
  p.dat[[k]] <- fit_sum[2,4]
  
}

dat.run <- B
dat.gp <- factor(phenoDataPAC12$Trimester)
dat.co1 <- factor(phenoDataPAC12$ArrayDateReceived)
p.dat <- list()
# dat.co2 <- factor(phenoDataPlacenta$Fetal.Sex, levels = c("F", "M"))

library(broom)
system.time(
  temp <- tibble(
    probe = rownames(B)[1:1000],
    fit = lapply(probe, function(x){
      beta.glm(B[x,], dat.gp, dat.co1)
    }),
    coef = lapply(fit, tidy),
    p = vapply(coef, function(x){dplyr::filter(x, term == "XSecond")$p.value}, numeric(1))
    # p = vapply(fit, function(x){
    #   suppressWarnings(summary(x)["XSecond", "Pr(>|t|)"])
    #   }, numeric(1))
  )
)


apply(B[1:3,], MARGIN = 1, function(x){
  fit <- beta.glm(x,dat.gp, dat.co1) #, dat.co2
  summary(fit)[2,4]

  p.dat_df <- data.frame(matrix(unlist(p.dat), nrow=length(p.dat), byrow=T),
                         stringsAsFactors=FALSE) %>% 
    `colnames<-`("unadj.p")
  
  rownames(p.dat_df) <- rownames(dat.run)
  p.dat_df$adj.P <- p.adjust(p.dat_df$unadj.p,"BH")
  
  # save files:
  # saveRDS(fit2, file = here("Method_DNAme/rds/glmBeta_fit_T1T2.rds"))
  saveRDS(p.dat_df, file = here("Method_DNAme/rds/glmBeta_pval_T1T2.rds"))
  
  return(p.dat_df=p.dat_df)
}


############## Use functions ############################################################

mbm <- microbenchmark(
  RUVm = RUVm_dmp(EpicRGsetPAC12 = EpicRGsetPAC12, M = BM_list$Mnorm),
  Lm = Lm_dmp(phenoDataPAC12 = phenoDataPAC12, M = BM_list$Mnorm),
  glmBeta = glmBeta_dmp(phenoDataPAC12 = phenoDataPAC12, B = BM_list$betaNorm)
)

saveRDS(mbm, file = here("Method_DNAme/rds/microbenchmark_DMPs_out.rds"))

library(ggplot2)
autoplot(mbm)

