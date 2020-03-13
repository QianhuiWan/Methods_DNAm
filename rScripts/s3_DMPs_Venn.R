
# load packages and data
library(grid)
library(VennDiagram)
library(eulerr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

B <- readRDS(file = here("Method_DNAme/rds/Bval.rds"))
EpicRGsetPAC12 <- read.metharray.exp(targets = phenoData, extended = TRUE, force=TRUE)

## calculate the delta beta values
table(phenoDataPAC12$Trimester)
keepT1 <- colnames(B) %in% rownames(phenoDataPAC12[phenoDataPAC12$Trimester=="First",])

# set up a new df 
B_T2_T1 <- B %>%  as.data.frame() %>% 
  magrittr::inset("meanT1", value=rowMeans(subset(., select = keepT1), na.rm = TRUE)) %>% 
  magrittr::inset("meanT2", value=rowMeans(subset(., select = !keepT1), na.rm = TRUE)) %>% 
  mutate(deltaBeta = meanT2-meanT1) %>% 
  `rownames<-`(rownames(B)) %>% 
  tibble::rownames_to_column(var="probe")

table(abs(B_T2_T1$deltaBeta) > 0.2)

# load DMP results and add delta beta to the DMPs result columns
limma <- readRDS(file = here("Method_DNAme/rds/Lm_dmp_res.rds")) %>% 
  tibble::rownames_to_column(var="probe") %>% 
  left_join(B_T2_T1[, c("probe", "deltaBeta")], by="probe")

RUVm_rinv <- readRDS(file = here("Method_DNAme/rds/RUVm_dmp_res.rds")) %>% 
  tibble::rownames_to_column(var="probe") %>% 
  left_join(B_T2_T1[, c("probe", "deltaBeta")], by="probe")

glmBeta_pval <- readRDS(file = here("Method_DNAme/rds/glmBeta_pval.rds")) %>% 
  left_join(B_T2_T1[, c("probe", "deltaBeta")], by="probe")
glmBeta_pval$p.adj <- p.adjust(glmBeta_pval$p, method = "fdr")

### compare methods <0.05:
# lmFit
## 24519
DMPs_lm_sig <- limma[limma$adj.P.Val < 0.05,]$probe

# RUVm (RUVrinv)
## 51284
DMPs_RUVm_rinv_sig <- RUVm_rinv[RUVm_rinv$p.BH_X1.Second<0.05, ]$probe

# Beta regression
## 51670
DMPs_glmBeta_sig <- glmBeta_pval[glmBeta_pval$p.adj<0.05, ]$probe

# intersect
intersect(DMPs_lm_sig,DMPs_RUVm_rinv_sig) %>% length() # 24508, 1-2
intersect(DMPs_RUVm_rinv_sig,DMPs_glmBeta_sig) %>% length() # 36927, 2-3
intersect(DMPs_lm_sig,DMPs_glmBeta_sig) %>% length() # 23653, 1-3

intersect(intersect(DMPs_lm_sig,DMPs_RUVm_rinv_sig), DMPs_glmBeta_sig)%>% length() # 23645


# Venn diagram

# input <- c(
#   lm=3,
#   RUVm=51284,
#   glmBeta=51670,
#   "lm&RUVm"=24508-23645,
#   "lm&glmBeta"=23653-23645,
#   "RUVm&glmBeta"=36927-23645,
#   "lm&RUVm&glmBeta"= 23645)
# 
# v <- euler(combinations = input,  shape = "ellipse") #shape="circle" ellipse
# 
# Venn_DMP_0.05 <- plot(v, 
#      # quantities=T,
#      quantities=list(col="black", cex=1),
#      # lty = 1:2,
#      # fills = list(fill = c("red", "steelblue4", "orange"), alpha = 0.8),
#      fills = list(fill = c("pink", "lightblue", "lightgoldenrodyellow"), alpha = 1),
#      # labels=F,
#      # labels = list(col = c("red", "steelblue4", "orange"), font=2, fontsize=20),
#      legend = list(labels = c("lm", "RUVm", "glmBeta"), 
#                    col= "black", alpha = 1, cex=2, side="right")) %>% as_ggplot()

Venn_DMP_0.05 <- draw.triple.venn(area1 = 24519,
                 area2 = 51284,
                 area3 = 51670,
                 n12 = 24508,
                 n23 = 36927,
                 n13 = 23653,
                 n123 = 23645,
                 category = c("lm", "RUVm", "glmBeta"),
                 fill = c("pink", "lightblue", "lightgoldenrodyellow"),
                 cex = 1.5, cat.cex = 1, 
                 cat.fontface = 2,
                 cat.dist = 0.04) 

Venn_DMP_0.05 <- grid.arrange(gTree(children=Venn_DMP_0.05), 
                              top=textGrob("FDR < 0.05",
                                           gp=gpar(fontface="bold",cex=1))) %>% 
  as_ggplot()

###########################################################################################
# use delta beta cut off to filter DMPs and draw the venn diagram again

### compare methods <0.05:
# lmFit
## 12725
DMPs_lm_2cf <- limma[limma$adj.P.Val < 0.05 & abs(limma$deltaBeta)>0.2,]$probe

# RUVm (RUVrinv)
## 19314
DMPs_RUVm_rinv_2cf <- RUVm_rinv[RUVm_rinv$p.BH_X1.Second<0.05 & abs(RUVm_rinv$deltaBeta)>0.2, ]$probe

# Beta regression
## 16182
DMPs_glmBeta_2cf <- glmBeta_pval[glmBeta_pval$p.adj<0.05& abs(glmBeta_pval$deltaBeta)>0.2, ]$probe

# intersect
intersect(DMPs_lm_2cf,DMPs_RUVm_rinv_2cf) %>% length() # 12720, 1-2
intersect(DMPs_RUVm_rinv_2cf,DMPs_glmBeta_2cf) %>% length() # 15745, 2-3
intersect(DMPs_lm_2cf,DMPs_glmBeta_2cf) %>% length() # 12636, 1-3

intersect(intersect(DMPs_lm_2cf,DMPs_RUVm_rinv_2cf), DMPs_glmBeta_2cf)%>% length() # 12631


# Venn diagram

Venn_DMP_2cf <- draw.triple.venn(
  area1 = 12725,
  area2 = 19314,
  area3 = 16182,
  n12 = 12720,
  n23 = 15745,
  n13 = 12636,
  n123 = 12631,
  category = c("lm", "RUVm", "glmBeta"),
  fill = c("pink", "lightblue", "lightgoldenrodyellow"),
  cex = 1.5, cat.cex = 1, 
  cat.fontface = 2,
  cat.dist = 0.04)

Venn_DMP_2cf <- grid.arrange(gTree(children=Venn_DMP_2cf), 
                             top=textGrob("FDR < 0.05 & ∆β ≥ 0.2",
                                          gp=gpar(fontface="bold",cex=1))) %>% 
  as_ggplot()




VennPlots_fig <- ggarrange(Venn_DMP_0.05, Venn_DMP_2cf,
         labels = c("A", "B"),
         font.label = list(size = 12, face = "plain", color ="black"), 
         nrow=1, ncol = 2)

jpeg(file = here("Method_DNAme/figures/VennPlots_DMPs.jpeg"),
     width = 16, height = 8, units = "cm", res = 320)
VennPlots_fig
dev.off()
