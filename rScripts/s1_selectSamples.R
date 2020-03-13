
library(here)
library(tidyverse)
library(magrittr)
library(minfi)

# read in the phenodata for all PAC placenta samples

baseDir <- file.path(here("EPICdata/IDAT_Placenta_PAC"))

phenoData <- read.metharray.sheet(baseDir)

phenoData_fil <- phenoData %>% 
  mutate(SampleName = paste0(Slide, "_", Array)) %>% 
  filter(Smoking.Status=="N") %>% 
  filter(Gestational.Age%in% c(6:10, 13:17)) %>% 
  filter(Ethnicity=="Caucasian") %>% 
  filter(Outliers%in%"Standard") %>% 
  filter(!(ArrayDateReceived=="2017-08-01T00:00:00Z")) %>% 
  filter(Sample.ID%in%c("PAC0010", "PAC0017", "PAC0026", "PAC0031", 
                        "PAC0032", "PAC0055", "PAC0064", "PAC0110",
                        "PAC0121", "PAC0148", "PAC0193", "PAC0198")) %>%  # "PAC0193" should be Male
  select(-Basename)

write_csv(phenoData_fil, path = here("EPICdata/IDAT_Placenta_PAC_sub12/placentaPhenoData_sub12_PAC.csv"))

# write 12 sample names to a txt file, so we can move the correspounding file into subPAC folder using bash scripts
phenoData_fil %>% 
  select(SampleName) %>% 
  write.table(file = here("PrepareData/moveIDATfiles/placenta_subPAC12_SampleName.txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)


# save files for samples
phenoData_fil %>% 
  mutate(Batch=str_replace_all(ArrayDateReceived, 
               c("2017-03-23T00:00:00Z"="1", "2019-02-01T00:00:00Z"="2"))) %>% 
  select(Tissue, Batch,Trimester, `Fetal sex`=Fetal.Sex) %>% 
  write_csv(path = here("Method_DNAme/files/Table3_methodPaper.csv"))
  

