
library(here)
library(tidyverse)
library(magrittr)
library(minfi)
library(ENmix)
library(ChAMP)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# read in the RGset ###########################################################
baseDir <- file.path(here("EPICdata/GSE126672/GSE126672_RAW"))

phenoData <- read.metharray.sheet(baseDir) %>% 
  dplyr::mutate(SampleName = paste0(Slide, "_", Array),`GEO accession`=GEO_accession, `Cell line`= Sample) %>% 
  dplyr::filter(str_detect(Sample, pattern = "WI-38"))

EpicRGsetGEO <- read.metharray.exp(targets = phenoData, extended = TRUE)
# phenoDataGEO <- pData(EpicRGsetPAC12) # chain data & phenotypes together

# preprocess ######################################################################

## QC ####################################################################
### PDF
qcReport(EpicRGsetGEO,
         sampNames = pData(EpicRGsetGEO)$`GEO accession`,
         sampGroups = pData(EpicRGsetGEO)$`Cell line`,
         pdf = here("Method_DNAme/files/", paste0("qcReport_GEO_4.pdf")))

### JPEG
out <- EpicRGsetGEO %>% preprocessRaw() %>% minfiQC

jpeg(file = here("Method_DNAme/figures/", paste0("plotQC_GEO.jpeg")),
     width = 20, height = 12, units = "cm", res = 300)
plotQC(out$qc)
dev.off()

### getSex and update phenoDataGEO
estSex <- EpicRGsetGEO %>% preprocessRaw() %>% mapToGenome() %>% getSex()
phenoDataGEO <- EpicRGsetGEO %>% preprocessRaw() %>% mapToGenome() %>% 
  addSex(sex = estSex) %>% 
  pData() # chain data & phenotypes together

jpeg(file = here("Method_DNAme/figures/", paste0("plotSex_GEO.jpeg")),
     width = 20, height = 12, units = "cm", res = 300)
EpicRGsetGEO %>% preprocessRaw() %>% mapToGenome() %>% 
  addSex(sex = estSex) %>% plotSex()
dev.off()

## filtering ###############################################################
crossreactiveProbes <- read_csv(here("PrepareData/DownloadedFilesForEPICpipeline","13059_2016_1066_MOESM1_ESM-cross-reactive.csv"))
colnames(crossreactiveProbes)[1] <- "ProbeNames"

islandHMM <- readRDS(file = here("PrepareData/DownloadedFilesForEPICpipeline", "HMM_CpGisland_hg19.rds"))

methylset_raw <- EpicRGsetGEO %>% preprocessRaw() # 865859

## drop P
detP <- minfi::detectionP(EpicRGsetGEO)
Failed <- detP > 0.01
methylset_raw_P <- methylset_raw[rowSums(Failed)==0,] # 860982

## drop probes with <3 beads

BeadCount <- getNBeads(EpicRGsetGEO)
RemainProbe <- rowSums(BeadCount<3) < 0.05 * (ncol(BeadCount))
methylset_rawBead <- EpicRGsetGEO[RemainProbe,] %>% preprocessRaw() # bead number filter

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
# rmXY <- !(featureNames(methylset_raw_PBeadXsnp) %in% ann850k$Name[ann850k$chr %in% c("chrX", "chrY")])
# methylset_raw_PBeadXsnpXY <- methylset_raw_PBeadXsnp[rmXY, ] # 606503

methylset_raw_PBeadXsnpXY <- methylset_raw_PBeadXsnp

# bg correction
filteredCpGs <- setdiff(featureNames(methylset_raw), featureNames(methylset_raw_PBeadXsnpXY))
methylset_bgCorr <- preprocessENmix(EpicRGsetGEO, exCpG = filteredCpGs, nCores = 3) 

## Nomalisation ####################################################################

beta_bgCorr <- minfi::getBeta(methylset_bgCorr)

betaNorm <- champ.norm(beta=beta_bgCorr,
                           # resultsDir="/home/qianhui/DNAme/Process_PAC_workflow/BMIQplots/",
                           method="BMIQ",
                           plotBMIQ=FALSE,
                           arraytype="EPIC",
                           cores=3)
### calculate M values
M <- log2(betaNorm/(1 - betaNorm))
Mnorm <- M[apply(M, 1, function(x) all(is.finite(x))), ] 

### Beta density plot
jpeg(file = here("Method_DNAme/figures/", 
                 paste0("density_normB_GEO.jpeg")),
     width = 20, height = 12, units = "cm", res = 300)
par(cex=0.8)
densityPlot(betaNorm,
            sampGroups = phenoDataGEO$`Cell line`, legend=FALSE,
            xlab = "Normalised beta values")
dev.off()


### PCA
#### PCA plot for normalised M values
PCA_M <- FactoMineR::PCA(base::t(Mnorm), graph = FALSE)

#### relabel PCA_PAC_M plot using ggplot:
PCA_M_plot.df <- data.frame(x = PCA_M$ind$coord[,1], y = PCA_M$ind$coord[,2],
                            Sex=phenoDataGEO$predictedSex,
                            Cell=phenoDataGEO$`Cell line`)
PCA_M_plot <- PCA_M_plot.df %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(aes(color= Cell), alpha = 0.9, size = 2)+ #color = factor(Batch),shape
  labs(title = paste0("PCA of normalised data"), x="PC1 (52.87%)", y="PC2 (28.46%)")

#### save plot
print(PCA_M_plot)
ggsave(here("Method_DNAme/figures/", paste0("PCA_M_plot_GEO.jpeg")),
       width = 20, height = 12, units = "cm", dpi = 300)



savelist <- list(betaNorm=betaNorm,
                 Mnorm=Mnorm,
                 phenoDataGEO=phenoDataGEO)

saveRDS(savelist, file = here("Method_DNAme/rds/preprocessed_GEO4.rds"))

