rm(list=ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(openxlsx)
  library(dplyr)
  library(tidyverse)
  library(enrichR)
  library(data.table)
})
setwd("/dsigdb")
input_file <- "Prioritized_eGenes.xlsx"
drug_gene <- read.xlsx(input_file)
if("Gene" %in% colnames(drug_gene)){symbol <- unique(drug_gene$Gene)
}else if("gene" %in% colnames(drug_gene)){symbol <- unique(drug_gene$gene)
}else{stop("Gene column not found")}
symbol <- symbol[!is.na(symbol)]
symbol <- unique(symbol)
cat("Input genes:",length(symbol),"\n")
############################
## DSigDB enrichment
############################
dbs <- listEnrichrDbs()
dbs <- dbs$libraryName[grepl("DSigDB",dbs$libraryName)]
enrich_result <- enrichr(symbol,dbs)
result <- enrich_result[[1]]
drug_result <- result %>% filter(P.value <0.05) %>% arrange(desc(Combined.Score))
drug_result$Target_Genes <- sapply(drug_result$Genes,
  function(x){length(unlist(strsplit(x,";")))})
drug_result$FDR <- p.adjust(drug_result$P.value,method="BH")
drug_unique <- drug_result %>% distinct(Term,.keep_all=TRUE)
top_drugs <- drug_unique %>% filter(FDR<0.05) %>% arrange(desc(Combined.Score))
write.xlsx(drug_result,"DSigDB_all_drug_results.xlsx",rowNames=FALSE)
write.xlsx(drug_unique,"DSigDB_unique_drug_results.xlsx",rowNames=FALSE)
write.xlsx(top_drugs,"DSigDB_candidate_drugs.xlsx",rowNames=FALSE)

