rm(list=ls())
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
  library(vcfR)
  library(openxlsx)
  library(stringr)
})
setwd("/DICE")
############################
## Parameters
############################
dice_files <- c(
  "B_CELL"="B_CELL_NAIVE.vcf",
  "CD4_T"="CD4_NAIVE.vcf",
  "CD8_T"="CD8_NAIVE.vcf",
  "MONOC"="MONOC.vcf",
  "MONONC"="MONONC.vcf",
  "NK"="NK.vcf"
)
candidate_file <- "coloc_PP_H4_gt_0.8.csv"
outcome_file <- "outcome_df.tsv"
############################
## Read candidate genes
############################
candidate <- fread(candidate_file)
candidate_genes <- unique(candidate$Gene)
############################
## Read AD outcome
############################
outcome_all <- fread(outcome_file)
############################
## DICE VCF processing function
############################
read_dice_vcf <- function(vcf_file,cell){
  vcf <- read.vcfR(vcf_file)
  fix <- as.data.frame(getFIX(vcf))
  info <- as.data.frame(getINFO(vcf))
eqtl_df <- bind_cols(
  fix,tibble(info=info[[1]]))%>%
  separate(info,into=c("Gene","GeneSymbol","Pvalue","Beta"),
  sep=";",extra="drop",
  fill="right")%>%
  mutate(Gene=sub("^Gene=","",Gene),
  GeneSymbol=sub("^GeneSymbol=","",GeneSymbol),
  Pvalue=as.numeric(sub("^Pvalue=","",Pvalue)),
  Beta=as.numeric(sub("^Beta=","",Beta)),
  Cell=cell)return(eqtl_df)}
############################
## MR replication function
############################
run_dice_mr <- function(eqtl_df,gene,cell){
dat <- eqtl_df %>% filter(GeneSymbol==gene)
if(nrow(dat)==0){return(NULL)}
exp_dat <- dat %>%
    transmute(
      SNP=ID,
      chr=CHROM,
      pos=POS,
      effect_allele=ALT, 
      other_allele=REF,
      beta=Beta,
      pval=Pvalue,exposure=paste0(gene,"_",cell))
exp_dat <- exp_dat %>%
mutate(z=qnorm(pval/2,lower.tail=FALSE),se=abs(beta)/z)%>% select(-z)
exp_dat <- exp_dat %>% filter(pval<5e-5)
if(nrow(exp_dat)==0){return(NULL)}
exp_dat <- format_data(
    exp_dat,
    type="exposure",
    snp_col="SNP",
    beta_col="beta",
    se_col="se",
    effect_allele_col="effect_allele",
    other_allele_col="other_allele",
    pval_col="pval",
    chr_col="chr",
    pos_col="pos",
    phenotype_col="exposure")
exp_dat <- clump_data(
    exp_dat,
    clump_kb=10000,
    clump_r2=0.001,
    clump_p1=1,clump_p2=1,pop="EUR")
if(nrow(exp_dat)==0){return(NULL)}
exp_dat <- exp_dat %>% mutate(F_stat= (beta.exposure^2)/(se.exposure^2))
if(min(exp_dat$F_stat,na.rm=TRUE)<=10){return(NULL)}
outcome <- outcome_all %>% filter(SNP %in% exp_dat$SNP)
if(nrow(outcome)==0){return(NULL)}
outcome_dat <- format_data(
    outcome,
    type="outcome",
    snp_col="SNP",
    beta_col="beta_outcome",
    se_col="se_outcome",
    effect_allele_col="ALT",
    other_allele_col="REF",
    pval_col="pval_outcome",
    phenotype_col="AD")
harm <- harmonise_data(exp_dat,outcome_dat,action=2)
harm <- harm %>% filter(mr_keep==TRUE)
if(nrow(harm)==0){}
n_snp <- nrow(harm)
if(n_snp==1){method="mr_wald_ratio"}else{method="mr_ivw"}    
mr_res <- mr(harm,method_list=method)
if(n_snp>=3){wm <- mr(harm,method_list="mr_weighted_median")
mr_res <- bind_rows(mr_res,wm)}
mr_res$Gene <- gene
mr_res$Cell <- cell
mr_res$N_SNP <- n_snp
mr_res$min_F <- min(harm$F_stat,na.rm=TRUE)
return(mr_res)}
############################
## Run all DICE datasets
############################
all_results <- list()
k <- 1
for(cell in names(dice_files)){cat("Processing:",cell,"\n")
eqtl <- read_dice_vcf(dice_files[cell],cell)
eqtl <- eqtl %>% filter(GeneSymbol %in% candidate_genes)
  for(gene in candidate_genes){
    res <- run_dice_mr(eqtl,gene,cell)
    if(!is.null(res)){all_results[[k]] <- res
      k <- k+1}}}
dice_results <- bind_rows(all_results)
dice_results$pval_bonf <- p.adjust(dice_results$pval,method="bonferroni")
dice_results <- dice_results %>%
    mutate(OR=exp(b),
    lower_CI=exp(b-1.96*se),
    upper_CI=exp(b+1.96*se),
    OR_95CI=sprintf("%.3f (%.3f-%.3f)",OR,lower_CI,upper_CI ))
############################
## Significant replication
############################
dice_sig <- dice_results %>% filter(pval_bonf<0.05)
############################
## Save
############################
write.xlsx(dice_results,"DICE_replication_all_results.xlsx",rowNames=FALSE)
write.xlsx(dice_sig,"DICE_replication_significant.xlsx",rowNames=FALSE)
cat("Total results:",nrow(dice_results),"\n")
cat("Significant:",nrow(dice_sig),"\n")