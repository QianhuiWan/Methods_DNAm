
library(here)
library(tidyverse)
library(magrittr)
library(minfi)
library(tidyverse); library(magrittr); library(GEOquery); library(here); library(minfi)

library(doParallel); library(minfi); library(IlluminaHumanMethylation27kmanifest)
library(ChAMP); library(wateRmelon); library(ggpubr); library(impute); library(WGCNA)
registerDoParallel(cores = 3)
library(RColorBrewer)
pal <- brewer.pal(8, "Dark2")

# get the truth set from GEO database
## get the raw data for GSE124565:
# getGEOSuppFiles("GSE126672", baseDir  = here("EPICdata"))

GSE126672 <- getGEO("GSE126672", destdir = here("EPICdata"))

## get meta data for this dataset 
GSE126672.meta <- pData(phenoData(GSE126672$GSE126672_series_matrix.txt.gz)) %>% 
  mutate(GEO_accession=geo_accession) %>% 
  mutate(Sample=`tissue:ch1`) %>% 
  mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
  mutate(Sample_nameToRowname=Sample_name) %>% 
  inset("Disease", value=rep('Cell line', nrow(.))) %>% 
  inset('Gestation', value=rep('NA', nrow(.))) %>% 
  inset('Trimester', value=rep('NA', nrow(.))) %>% 
  mutate(Additional_Info=source_name_ch1) %>% 
  mutate(Additional_Info_2=`description`) %>% 
  inset("Fetal_Sex", value= rep('NA', nrow(.))) %>% 
  inset('ArrayType', value=rep('EPIC', nrow(.))) %>% 
  inset('Study', value=rep('GSE126672', nrow(.))) %>%
  mutate(Array= str_extract(Sample_name, pattern = "R..C..")) %>% 
  mutate(Slide=  str_extract(Sample_name, pattern= "_(.*?)_")) %>% 
  mutate(Slide=  str_replace_all(Slide, '_', '')) %>% 
  dplyr::select(GEO_accession:Slide) %>% 
  column_to_rownames(var='Sample_nameToRowname')

write_csv(GSE126672.meta, path = here("EPICdata/GSE126672/GSE126672_RAW/GSE126672.meta.csv"))

# read in the phenodata for all PAC placenta samples
baseDir <- file.path(here("EPICdata/GSE126672/GSE126672_RAW"))

phenoData <- read.metharray.sheet(baseDir)

phenoData %>% 
  dplyr::mutate(SampleName = paste0(Slide, "_", Array)) %>% 
  dplyr::filter(str_detect(Sample, pattern = "WI-38")) %>% 
  dplyr::select(-Basename) %>% 
  dplyr::select(`GEO accession`=GEO_accession, `Cell line`= Sample) %>% 
  write_csv(path = here("Method_DNAme/files/Table_methodPaper_truthDataSet.csv"))
  