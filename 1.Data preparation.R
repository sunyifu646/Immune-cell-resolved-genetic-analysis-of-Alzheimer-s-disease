#############################################################
# OneK1K sc-eQTL data preparation
# Input:
# OneK1K cell-type-specific cis-eQTL summary statistics
#
# Output:
# Final independent IV SNPs
#############################################################

rm(list = ls())

library(data.table)
library(dplyr)
library(ieugwasr)
library(stringr)

setwd("~/OneK1K_snp/")
cell_files <- c(
  
  "B_IN"="bin_eqtl_table.tsv.gz",
  
  "B_MEM"="bmem_eqtl_table.tsv.gz",
  
  "Plasma"="plasma_eqtl_table.tsv.gz",
  
  "CD4_NC"="cd4nc_eqtl_table.tsv.gz",
  
  "CD4_ET"="cd4et_eqtl_table.tsv.gz",
  
  "CD4_SOX4"="cd4sox4_eqtl_table.tsv.gz",
  
  "CD8_NC"="cd8nc_eqtl_table.tsv.gz",
  
  "CD8_ET"="cd8et_eqtl_table.tsv.gz",
  
  "CD8_S100B"="cd8s100b_eqtl_table.tsv.gz",
  
  "NK"="nk_eqtl_table.tsv.gz",
  
  "NKR"="nkr_eqtl_table.tsv.gz",
  
  "Mono_C"="monoc_eqtl_table.tsv.gz",
  
  "Mono_NC"="mononc_eqtl_table.tsv.gz",
  
  "DC"="dc_eqtl_table.tsv.gz"
  
)
process_eqtl <- function(file, cell){
  
  
  cat("\nProcessing:",cell,"\n")
  
  
  eqtl <- fread(
    paste0(
      "data/OneK1K/",
      file
    )
  )
  
  
  
  cat(
    "Original SNP:",
    nrow(eqtl),
    "\n"
  )
  
  
  
################################################
  # genome-wide significant cis-eQTL
################################################
  
  
  eqtl <- eqtl %>%
    
    filter(
      pval < 5e-8
    )
  
  
  
  cat(
    "p<5e-8:",
    nrow(eqtl),
    "\n"
  )
  
  
  
################################################
  # standardize column names
################################################
  
  
  eqtl <- eqtl %>%
    
    rename(
      
      SNP = variant_id,
      
      Gene = gene_id,
      
      Beta = beta,
      
      SE = se,
      
      P = pval,
      
      EA = effect_allele,
      
      OA = other_allele
      
    )
  
  
  
  eqtl$Cell <- cell
  
  
  
  return(eqtl)
  
  
}
all_eqtl <- list()



for(i in names(cell_files)){
  
  
  all_eqtl[[i]] <-
    
    process_eqtl(
      
      cell_files[i],
      
      i
      
    )
  
  
}


all_eqtl <- rbindlist(
  all_eqtl
)


dim(all_eqtl)

################################################
# LD clumping
################################################
clump_input <- data.frame(
  
  rsid = all_eqtl$SNP,
  
  pval = all_eqtl$P,
  
  id = paste0(
    all_eqtl$Gene,
    "_",
    all_eqtl$Cell
  )
  
)

clumped <- ld_clump(
  
  dat = clump_input,
  
  
  clump_kb = 10000,
  
  
  clump_r2 = 0.001,
  
  
  pop="EUR"
  
  
)



length(
  unique(clumped$rsid)
)
##
eqtl_independent <- all_eqtl %>%
  
  filter(
    
    SNP %in%
      
      clumped$rsid
    
  )
eqtl_independent <- eqtl_independent %>%
  
  
  mutate(
    
    F_stat=
      
      (Beta^2)/(SE^2)
    
  )
eqtl_final <- eqtl_independent %>%
  
  filter(
    
    F_stat >=10
    
  )
######
fwrite(
  
  eqtl_final,
  
  "result/OneK1K_final_SNP.txt"
  
)
######
cell_summary <- eqtl_final %>%
  
  group_by(Cell) %>%
  
  summarise(
    
    SNP_number =
      n_distinct(SNP),
    
    Gene_number =
      n_distinct(Gene)
    
  )


write.csv(
  
  cell_summary,
  
  "result/OneK1K_cell_summary.csv",
  
  row.names=FALSE
  
)