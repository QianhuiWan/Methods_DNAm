
library(here)
library(tidyverse)
library(magrittr)
library(minfi)
library(ENmix)
library(ChAMP)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# read in the RGset ###########################################################
baseDir <- file.path(here("EPICdata/IDAT_Placenta_PAC_sub12"))

phenoData <- read.metharray.sheet(baseDir)

EpicRGsetPAC12 <- read.metharray.exp(targets = phenoData, extended = TRUE, force = TRUE)
# phenoDataPAC12 <- pData(EpicRGsetPAC12) # chain data & phenotypes together

# preprocess ######################################################################

## QC ####################################################################
### PDF
qcReport(EpicRGsetPAC12,
         sampNames = pData(EpicRGsetPAC12)$Sample.ID,
         sampGroups = pData(EpicRGsetPAC12)$Trimester,
         pdf = here("Method_DNAme/files/", paste0("qcReport_placenta_PAC_12.pdf")))

### JPEG
out <- EpicRGsetPAC12 %>% preprocessRaw() %>% minfiQC

jpeg(file = here("Method_DNAme/figures/", paste0("plotQC_PAC12_placenta.jpeg")),
     width = 20, height = 12, units = "cm", res = 300)
plotQC(out$qc)
dev.off()

### getSex and update phenoDataPAC12
estSex <- EpicRGsetPAC12 %>% preprocessRaw() %>% mapToGenome() %>% getSex()
phenoDataPAC12 <- EpicRGsetPAC12 %>% preprocessRaw() %>% mapToGenome() %>% 
  addSex(sex = estSex) %>% 
  pData() # chain data & phenotypes together

jpeg(file = here("Method_DNAme/figures/", paste0("plotSex_PAC12_placenta.jpeg")),
     width = 20, height = 12, units = "cm", res = 300)
EpicRGsetPAC12 %>% preprocessRaw() %>% mapToGenome() %>% 
  addSex(sex = estSex) %>% plotSex()
dev.off()

## filtering ###############################################################
crossreactiveProbes <- read_csv(here("PrepareData/DownloadedFilesForEPICpipeline","13059_2016_1066_MOESM1_ESM-cross-reactive.csv"))
colnames(crossreactiveProbes)[1] <- "ProbeNames"

islandHMM <- readRDS(file = here("PrepareData/DownloadedFilesForEPICpipeline", "HMM_CpGisland_hg19.rds"))

methylset_raw <- EpicRGsetPAC12 %>% preprocessRaw() # 865859

## drop P
detP <- minfi::detectionP(EpicRGsetPAC12)
Failed <- detP > 0.01
methylset_raw_P <- methylset_raw[rowSums(Failed)==0,] # 860982

## drop probes with <3 beads

BeadCount <- getNBeads(EpicRGsetPAC12)
RemainProbe <- rowSums(BeadCount<3) < 0.05 * (ncol(BeadCount))
methylset_rawBead <- EpicRGsetPAC12[RemainProbe,] %>% preprocessRaw() # bead number filter

methylset_raw_PBead <- methylset_raw_P[rownames(methylset_raw_P)%in%rownames(methylset_rawBead),] # 820854

## drop cross-reactive probes
keep <- !(featureNames(methylset_raw_PBead) %in% crossreactiveProbes$ProbeNames)
methylset_raw_PBeadX <- methylset_raw_PBead[keep, ] # 779489

## drop SNPs
GRaSet_raw_PBeadXsnp <- methylset_raw_PBeadX %>% 
  ratioConvert() %>% mapToGenome() %>% 
  addSnpInfo() %>% 
  dropLociWithSnps(snps = c("CpG", "SBE", "Probe"), maf = 0)

methylset_raw_PBeadXsnp <-
  methylset_raw_PBeadX[rownames(methylset_raw_PBeadX) %in% rownames(GRaSet_raw_PBeadXsnp), ] # 621684

# saveRDS(methPAC125_PBeadXsnp, file = here())

## drop probes on X, Y chr
rmXY <- !(featureNames(methylset_raw_PBeadXsnp) %in% ann850k$Name[ann850k$chr %in% c("chrX", "chrY")])
methylset_raw_PBeadXsnpXY <- methylset_raw_PBeadXsnp[rmXY, ] # 606503
# saveRDS(methPAC125_PBeadXsnpXY, file = here())

# bg correction
filteredCpGs <- setdiff(featureNames(methylset_raw), featureNames(methylset_raw_PBeadXsnpXY))
methylset_bgCorr <- preprocessENmix(EpicRGsetPAC12, exCpG = filteredCpGs, nCores = 3) 

## Nomalisation ####################################################################

beta_bgCorr <- minfi::getBeta(methylset_bgCorr)

betaNorm <- champ.norm(beta=beta_bgCorr,
                           # resultsDir="/home/qianhui/DNAme/Process_PAC_workflow/BMIQplots/",
                           method="BMIQ",
                           plotBMIQ=FALSE,
                           arraytype="EPIC",
                           cores=2)
### calculate M values
M <- log2(betaNorm/(1 - betaNorm))
Mnorm <- M[apply(M, 1, function(x) all(is.finite(x))), ] 

### Beta density plot
jpeg(file = here("Method_DNAme/figures/", 
                 paste0("density_normB_placentaPAC12.jpeg")),
     width = 20, height = 12, units = "cm", res = 300)
par(cex=0.8)
densityPlot(betaNorm,
            sampGroups = phenoDataPAC12$Trimester, legend=TRUE,
            xlab = "Normalised beta values")
dev.off()

### SVD
pd <- phenoDataPAC12[, c("ArrayDateReceived", "Gestational.Age", "Trimester", 
                         "Oxygen", "Fetal.Sex", "Maternal.Age", 
                         "BMI", "Smoking.Status", "Ethnicity")]

## chage to right class for each column
cols = c("Gestational.Age","Maternal.Age", "BMI")
pd[,cols] <-  lapply(pd[, cols], function(x){as.numeric(x)})
## champ.SVD
champ.SVD(beta = betaNorm, pd = pd, 
          resultsDir = here("Method_DNAme/figures/BeforeBatchCorrection"))


## Batch correction #######################################################################

combatAdjusted_Beta <- champ.runCombat(beta=betaNorm,
                                           pd=pd,
                                           variablename="Oxygen",
                                           batchname=c("ArrayDateReceived"),
                                           logitTrans=TRUE)

## calculate conbat adjusted M values
combat_Mnorm <- log2(combatAdjusted_Beta/(1 - combatAdjusted_Beta))
combatAdjusted_Mnorm <- combat_Mnorm[apply(combat_Mnorm, 1, function(x) all(is.finite(x))), ] 

### PCA
#### PCA plot for normalised M values
PCA_M <- FactoMineR::PCA(base::t(combatAdjusted_Mnorm), graph = FALSE)

#### relabel PCA_PAC_M plot using ggplot:
PCA_M_plot.df <- data.frame(x = PCA_M$ind$coord[,1], y = PCA_M$ind$coord[,2],
                            fetalSex=phenoDataPAC12$Fetal.Sex,
                            gestationalAge=phenoDataPAC12$Gestational.Age,
                            Batch=phenoDataPAC12$ArrayDateReceived,
                            trimester=phenoDataPAC12$Trimester)
PCA_M_plot <- PCA_M_plot.df %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(aes( color= Batch), alpha = 0.9, size = 2)+ #color = factor(Batch),shape
  geom_text(aes(label=factor(gestationalAge)), hjust=0, vjust=0)+ #, col=trimester
  labs(title = paste0("PCA of normalised data"), x="PC1", y="PC2")

#### save plot
print(PCA_M_plot)
ggsave(here("Method_DNAme/figures/", paste0("PCA_M_plot_placentaPAC12.jpeg")),
       width = 20, height = 12, units = "cm", dpi = 300)

### SVD analysis after Combat batch correction
champ.SVD(beta = combatAdjusted_Beta, pd = pd, 
          resultsDir = here("Method_DNAme/figures/AfterBatchCorrection"))


savelist <- list(betaNorm=betaNorm,
                 Mnorm=Mnorm,
                 # Adjusted_Beta=combatAdjusted_Beta,
                 # Adjusted_Mnorm=combatAdjusted_Mnorm,
                 phenoDataPAC12=phenoDataPAC12)

saveRDS(savelist, file = here("Method_DNAme/rds/preprocessed_PAC12_placenta.rds"))

