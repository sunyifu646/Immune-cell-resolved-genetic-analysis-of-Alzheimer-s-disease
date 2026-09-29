rm(list=ls())
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(coloc)
  library(locuscomparer)
})
setwd("")
############################
## Parameters
############################
eqtl_file <- "OneK1K_final_SNP.txt"
mr_file <- "Sc_eQTL_MR_significant_results.csv"
outcome_file <- "outcome_df.tsv"
window <- 90000
N_eqtl <- 928
N_outcome <- 487511
N_cases <- 85934
s_value <- N_cases/N_outcome
############################
## Read data
############################
eqtl <- fread(eqtl_file)
mr_sig <- fread(mr_file)
outcome_df <- fread(outcome_file)
############################
## Check columns
############################
eqtl_cols <- c(
  "SNP",
  "Gene",
  "Cell",
  "Beta",
  "SE",
  "EA",
  "OA",
  "P"
)
if(!all(eqtl_cols %in% colnames(eqtl))){
  stop("Missing columns in eQTL file")
}
outcome_cols <- c(
  "SNP",
  "CHR",
  "POS",
  "beta_outcome",
  "se_outcome",
  "maf_outcome",
  "pval_outcome"
)
if(!all(outcome_cols %in% colnames(outcome_df))){
  stop("Missing columns in outcome file")
}
mr_cols <- c(
  "Gene",
  "Cell")
if(!all(mr_cols %in% colnames(mr_sig))){
  stop("Missing columns in MR result")
}
############################
## Select MR significant gene-cell pairs
############################
target_pairs <- mr_sig %>%
  select(
    Gene,
    Cell
  ) %>%
  distinct()
cat("Number of target gene-cell pairs:",nrow(target_pairs),
  "\n")
############################
## Output directories
############################
dir.create(
  "coloc_result",
  showWarnings=FALSE)
dir.create(
  "coloc_result/locuscompare",
  showWarnings=FALSE)
############################
## Generate coloc regions
############################
run_coloc <- function(gene,cell){
  cat("Running:",gene,cell,"\n")
  eqtl_region <- eqtl %>%
    filter(Gene==gene,Cell==cell)
  if(nrow(eqtl_region)==0){
    return(NULL)}
  lead <- eqtl_region %>%
    arrange(P) %>%
    slice(1)
  chr <- lead$CHR
  pos <- lead$POS
  start <- pos-window
  end <- pos+window
  eqtl_region <- eqtl_region %>%
    filter(CHR==chr,POS>=start,POS<=end)
  if(nrow(eqtl_region)<10){
   cat("Skip: insufficient eQTL SNPs\n")
   return(NULL)}
  eqtl_region <- eqtl_region %>%
    select(
      SNP,
      CHR,
      POS,
      Beta,
      SE,
      P,
      EA,
      OA
    ) %>%
  distinct(SNP,.keep_all=TRUE)
  outcome_region <- outcome_df %>%
    filter(CHR==chr,POS>=start,POS<=end) %>%
    distinct(SNP,.keep_all=TRUE)
  merged <- inner_join(
    eqtl_region,
    outcome_region,
    by="SNP")
  if(nrow(merged)<10){cat("Skip: insufficient shared SNPs\n")
    return(NULL)}
dataset1 <- list(
    snp=merged$SNP,
    beta=merged$Beta,
    varbeta=merged$SE^2,
    MAF=NA,
    N=N_eqtl,
    type="quant")
dataset2 <- list(
    snp=merged$SNP,
    beta=merged$beta_outcome,
    varbeta=merged$se_outcome^2,
    MAF=merged$maf_outcome,
    N=N_outcome,
    type="cc",
    s=s_value)
coloc_res <- try(
    coloc.abf(
      dataset1,
      dataset2),
    silent=TRUE)
if(inherits(coloc_res,"try-error")){
cat("coloc failed\n")
    return(NULL)}
  summary <- coloc_res$summary
  result <- data.frame(
    Gene=gene,
    Cell=cell,
    CHR=chr,
    START=start,
    END=end,
    NSNP=as.numeric(summary["nsnps"]),
    PP_H0=as.numeric(summary["PP.H0.abf"]),
    PP_H1=as.numeric(summary["PP.H1.abf"]),
    PP_H2=as.numeric(summary["PP.H2.abf"]),
    PP_H3=as.numeric(summary["PP.H3.abf"]),
    PP_H4=as.numeric(summary["PP.H4.abf"]))
prefix <- paste0(gene,
    "_",gsub("[ /]","_",cell))
  saveRDS(
    coloc_res,
    paste0(
      "coloc_result/",
      prefix,
      "_coloc.rds"))
  eqtl_locus <- merged %>%
    select(
      SNP,
      P)
  colnames(eqtl_locus)<-c(
    "rsid",
    "pval")
  gwas_locus <- merged %>%
    select(SNP,pval_outcome)
  colnames(gwas_locus)<-c(
    "rsid",
    "pval")
  eqtl_file_tmp <- paste0(
    "coloc_result/",
    prefix,
    "_eqtl.tsv")
  gwas_file_tmp <- paste0(
    "coloc_result/",
    prefix,
    "_gwas.tsv")
  write.table(
    eqtl_locus,
    eqtl_file_tmp,
    sep="\t",
    quote=FALSE,
    row.names=FALSE)
  write.table(
    gwas_locus,
    gwas_file_tmp,
    sep="\t",
    quote=FALSE,
    row.names=FALSE)
  png(paste0(
      "coloc_result/locuscompare/",
      prefix,
      "_locuscompare.png"),
    width=3000,
    height=2000,
    res=300)
  locuscompare(
    in_fn1=eqtl_file_tmp,
    in_fn2=gwas_file_tmp,
    marker_col1="rsid",
    pval_col1="pval",
    marker_col2="rsid",
    pval_col2="pval",
    title1=paste0(gene," ",cell),
    title2="AD GWAS")
  dev.off()
  return(result)}
############################
## Run coloc
############################
coloc_list <- list()
for(i in seq_len(nrow(target_pairs))){
  res <- run_coloc(
    target_pairs$Gene[i],
    target_pairs$Cell[i])
  if(!is.null(res)){
    coloc_list[[length(coloc_list)+1]] <- res}}
############################
## Combine results
############################
coloc_summary <- bind_rows(coloc_list)
coloc_strong <- coloc_summary %>% filter(PP_H4>0.8)
############################
## Save
############################
fwrite(coloc_summary,"coloc_result/coloc_all_results.csv")
fwrite(coloc_strong,"coloc_result/coloc_PP_H4_gt_0.8.csv")
cat("Total coloc:",nrow(coloc_summary),"\n")
cat("Strong PP.H4 >0.8:",nrow(coloc_strong),"\n")