###############################################################
# Single-cell eQTL MR analysis
#
# Input:
# 1. OneK1K_final_SNP.txt
# 2. ebi-a-GCST90027158.vcf.gz
#
# Output:
# Immune cell-specific eGene MR results
###############################################################
rm(list = ls())
library(data.table)
library(dplyr)
library(stringr)
library(vcfR)
library(TwoSampleMR)
###############################################################
# Working directory
###############################################################
setwd("~/AD_sc_eQTL_MR/")
dir.create(
  "result",
  showWarnings = FALSE)
###############################################################
# 1. Read OneK1K final instruments
###############################################################
iv_dat <- fread("data/OneK1K/OneK1K_final_SNP.txt")
head(iv_dat)
###############################################################
# 2. Read AD GWAS VCF
###############################################################
vcf <- vcfR::read.vcfR("data/AD_GWAS/ebi-a-GCST90027158.vcf.gz")
gt <- data.frame(vcf@gt)
fix <- data.frame(vcf@fix)
parts <- strsplit(gt$ebi.a.GCST90027158,":",fixed = TRUE)
table(lengths(parts))
keep <- lengths(parts)==5
parts_keep <- parts[keep]
fix_keep <- fix[keep,]
mat <- do.call(rbind,parts_keep)
frame <- as.data.frame(mat, stringsAsFactors = FALSE)
colnames(frame) <- c(
  "ES",
  "SE",
  "LP",
  "AF",
  "ID")
AD <- cbind(fix_keep, frame)
AD <- AD %>% mutate(
    ES = as.numeric(ES),
    SE = as.numeric(SE),
    LP = as.numeric(LP),
    AF = as.numeric(AF),
    pval = 10^(-LP))
all_snps <- unique(iv_dat$SNP)
AD <- AD %>% filter(ID %in% all_snps)
###############################################################
# Create outcome dataset
###############################################################
outcome_dat <- AD %>%
  transmute(
    NP = ID,
    chr = CHROM,
    pos = POS,
    effect_allele.outcome = ALT,
    other_allele.outcome = REF,
    beta.outcome = ES,
    se.outcome = SE,
    eaf.outcome = AF,
    pval.outcome = pval,
    samplesize.outcome = 487511)
head(outcome_dat)
###############################################################
# 3.Create exposure dataset
###############################################################
make_exposure <- function(dat){
  exp <- dat %>% 
      transmute(
      SNP = SNP,
      effect_allele.exposure = EA,
      other_allele.exposure = OA,
      beta.exposure = Beta,
      se.exposure = SE,
      eaf.exposure = NA,
      pval.exposure = P,
      samplesize.exposure = NA)
  exp$id.exposure <- paste0(
    unique(dat$Gene),
    "_",
    unique(dat$Cell))exp$exposure <- exp$id.exposure
  return(exp)
} 
###############################################################
# 4. MR function
###############################################################
run_MR <- function(exp_dat,outcome_dat){
  harm <- harmonise_data(exp_dat,outcome_dat,action=2)
  harm <- harm %>% filter(mr_keep == TRUE)
  if(nrow(harm)==0){return(NULL)}
  nsnp <- nrow(harm)
  if(nsnp==1){res <- mr(harm,method_list=c("mr_wald_ratio"))}
  if(nsnp>=2){res <- mr(harm,method_list=c("mr_ivw"))}
  if(nsnp>=3){wm <- mr(harm,method_list=c("mr_weighted_median"))res <- bind_rows(res,wm)
  }
  res$nsnp <- nsnp
  return(res)
}
###############################################################
# 5. Run MR 
###############################################################
gene_cell <- iv_dat %>%
  group_by(
    Gene, Cell
  ) %>%
  group_split()
length(gene_cell)
MR_list <- list()
count <- 1
for(dat in gene_cell){
  gene <- unique(dat$Gene)
  cell <- unique(dat$Cell)
  cat("Running:",gene,cell,"\n")
  exp_dat <- make_exposure(dat)
  res <- run_MR(exp_dat,outcome_dat)
if(!is.null(res)){res$Gene <- gene
  res$Cell <- cell
  MR_list[[count]] <- res
  count <- count + 1
}
}
MR_results <- rbindlist(MR_list,fill=TRUE)
###############################################################
# 6. Bonferroni correction
###############################################################
MR_results <- MR_results %>% mutate(P_Bonferroni = pval * 8733)
MR_sig <- MR_results %>% filter(P_Bonferroni <0.05)
###############################################################
# 7. Save results
###############################################################
fwrite(MR_results,"result/Sc_eQTL_MR_all_results.csv")
fwrite(MR_sig,"result/Sc_eQTL_MR_significant_results.csv")
