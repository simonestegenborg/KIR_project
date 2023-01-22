################################### VDJ table preperation #######################################################
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(scoper)

# filtered_contig_annotations
# Load data 
vdj_data <- read.csv(file = "/Volumes/T-cells-and-cancer/SRH group/Group members/Simone/data/filtered_contig_annotations.csv")
vdj_data <- read.csv(file = "/Users/simone/Desktop/MASTER/metadata/filtered_contig_annotations.csv")

# Create one row pr. barcode and included only certain parameters
#separate parameters dependent on TCR alpha or beta chain by adding _TRA and _TRB
vdj_data_wrangled <- vdj_data %>% 
  pivot_wider(id_cols = barcode, 
              names_from = chain, 
              values_from = c(v_gene, d_gene, j_gene, c_gene, fwr1, fwr1_nt, cdr1, cdr1_nt, fwr2, fwr2_nt, cdr2, cdr2_nt, fwr3, fwr3_nt, cdr3, cdr3_nt, fwr4, fwr4_nt, umis), 
              names_sep = "_") 

# Remove d_gene_TRA
vdj_data_wrangled <- vdj_data_wrangled %>% select(-d_gene_TRA)

# Remove doublets completely
vdj_data_wrangled <- vdj_data_wrangled[!grepl(",", vdj_data_wrangled$v_gene_TRB),] 
vdj_data_wrangled <- vdj_data_wrangled[!grepl(",", vdj_data_wrangled$v_gene_TRA),] 

## Divide tables into two according to TRB and TRA
#splitting the table by suffix and splitting the large list created
vdj_data_seperated <- map(set_names(c("TRB", "TRA")),~select(vdj_data_wrangled,ends_with(.x)))
mapply(assign, names(vdj_data_seperated), vdj_data_seperated, MoreArgs=list(envir = globalenv()))

#binding the barcode column to each table
TRA <- cbind(vdj_data_wrangled[,1, drop=FALSE], TRA)
TRB <- cbind(vdj_data_wrangled[,1, drop=FALSE], TRB)


## Find identical clones
TRA[] <- lapply(TRA, as.character)
TRB[] <- lapply(TRB, as.character)

TRA <- identicalClones(TRA %>% 
                         mutate(locus_TRA="TRA") %>%
                         filter(!is.null(cdr3_nt_TRA) & cdr3_nt_TRA != "NULL"),
                       method="nt",junction = "cdr3_nt_TRA", 
                       v_call = "v_gene_TRA", j_call = "j_gene_TRA", summarize_clones = F,
                       clone="TRA_clone", cell_id=NULL, locus="locus_TRA", verbose = T)


TRB <- identicalClones(TRB %>% 
                         mutate(locus_TRB="TRB") %>%
                         filter(!is.null(cdr3_nt_TRB) & cdr3_nt_TRB != "NULL"),
                       method="nt", junction = "cdr3_nt_TRB", 
                       v_call = "v_gene_TRB", j_call = "j_gene_TRB", summarize_clones = F,
                       clone="TRB_clone", cell_id=NULL, locus="locus_TRB", verbose = T)


# Merge the two tables together and create an clone_id with "clones_TRB _ clones_TRA"
merged_vdj_table <- full_join(TRA,TRB,by="barcode")
merged_vdj_table <- merged_vdj_table %>% mutate(clone_id= paste(TRB_clone,TRA_clone,sep="_"))

# Change "NULL"s to NA
merged_vdj_table  <- merged_vdj_table %>% mutate_all(na_if,"NULL")

# Save as csv file 
write.csv(merged_vdj_table,file = "/Volumes/T-cells-and-cancer/SRH group/Group members/Simone/results/vdj_table.csv")


