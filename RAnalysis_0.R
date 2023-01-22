### Masters Project ###
library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)
library(patchwork) 
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(biomaRt)  
library(scRepertoire)
library(ggalluvial)
library(ggvenn)
library(ggVennDiagram)
library(venn)
library(ggplot2)
library(ggpolypath)
library(gprofiler2)
library(ComplexHeatmap)


#################### Create seuratobject and metadata ##################

vdj_data <- read.csv(file = "/Volumes/T-cells-and-cancer/SRH group/Group members/Simone/data/metadata/vdj_table.csv")
#vdj_data <- read.csv(file = "/Users/simone/Desktop/MASTER/vdj_table.csv")

#Creating a seurat object from the three files in /home/projects/SRHgroup/projects/SingleCell_KIR/results/filtered_cellranger_results_for_R
## barcodes.tsv.gz
## features.tsv.gz
## matrix.mtx.gz

# Initialize the Seurat object
data_dir <- '/Volumes/T-cells-and-cancer/SRH group/Group members/Simone/data/filtered_cellrranger_results'
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
KIRdata <- Read10X(data.dir = data_dir)
KIR <- CreateSeuratObject(counts = KIRdata$`Gene Expression`, project = "KIR")
KIR[['ADT']] = CreateAssayObject(counts = KIRdata$`Antibody Capture`)

# Add VDJ to metadata:
KIR <- AddMetaData(KIR, metadata = vdj_data %>% column_to_rownames("barcode"))

# Demultiplexing 
sample.list <- c("S1", "S2", "S3", "S4","S5","S6","S7","S8","S9","S10")
#extract sample data from Antibody Capture matrix
sample <- as.matrix(GetAssayData(object = KIR[["ADT"]], slot = "counts"))[sample.list,]
#create assay object 
KIR[["HTO"]] <- CreateAssayObject(counts = sample)


# Normalize and demultiplex
KIR <- KIR %>% 
  NormalizeData(assay = "HTO", normalization.method = "CLR") %>%
  HTODemux(assay = "HTO", positive.quantile = 0.99)

# Demultiplexing visualizations
#display singlet, doublet, and empty droplet count
table(KIR$HTO_classification.global)
table(KIR$HTO_classification)


#group cells based on the max HTO signal
Idents(KIR) <- "HTO_maxID"
RidgePlot(KIR, assay = "HTO", features = c("S1", "S2", "S3", "S4","S5","S6","S7","S8","S9","S10"))

# HTO Heatmap
KIR$HTO_classification2 <- KIR$HTO_classification
KIR$HTO_classification2[grepl("_", KIR$HTO_classification2)] <- "Multiplet"

Idents(KIR) <- factor(KIR$HTO_classification2, 
                       levels = c("S1","S2","S3","S4","S5",
                                  "S6","S7","S8","S9","S10","Multiplet", "Negative"))

KIR %>% subset(HTO_classification2 != "Negative") %>%
  DoHeatmap(assay = "HTO", slot = "data", group.bar = T,
          features = rownames(KIR@assays$HTO@counts), size = 3)  + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 11, vjust = 0.5)) +  
  guides(color='none') 


#Violingplot that compare number of UMIs for singlets, doublets and negative cells
Idents(KIR) <- "HTO_classification.global"
VlnPlot(KIR, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


# Create metadata stating whether each sample has KIR staining and buffy coat
KIR$KIR_staining <- "F"
KIR$buffy_coat   <- "BC351"

KIR$KIR_staining[KIR$hash.ID %in% c("S4","S5","S8","S9","S10")] <- "T"
KIR$buffy_coat[KIR$hash.ID %in% c("S6","S7","S8","S9","S10")] <- "BC357"

table(KIR$KIR_staining)
table(KIR$buffy_coat)

# If necessary, take only one buffy_coat here:
# KIR <- subset(KIR, buffy_coat == "BC3517")

# Label cells based on barcode
barcode.list <- c("BC1","BC2","BC3","BC4","BC5","BC6","BC7","BC9","BC10","BC11","BC12","BC13","BC14","BC15","BC16")

# Extract barcode matrix 
BC_raw <- as.data.frame(as.matrix(GetAssayData(object = KIR@assays$ADT, slot = "counts")))[barcode.list,]
rowSums(BC_raw)

# Exclude barcodes with low count to avoid zero-clusters and exclude GEMS with zero barcodes based on rowSums(BC_raw)
BC <- BC_raw[-c(8,9,10,11,12,14,15),]
BC <- BC[, colSums(BC) > 0]

# check barcode sums
rowSums(BC)

# Remove samples without those barcodes from seurat object to enable demultiplexing
demux_data <- subset(KIR, cells = colnames(BC))

# Add BC information to data object 
demux_data[['Barcodes']] <- CreateAssayObject(counts = BC)

# Demulitplex barcode information 
demux_data <- demux_data %>% 
  NormalizeData(assay = "Barcodes", normalization.method = "CLR") %>%
  HTODemux(assay = "Barcodes", positive.quantile = 0.99, seed = 46)

# Inspect data
table(demux_data$Barcodes_classification.global)
table(demux_data$Barcodes_classification)

# Add information as metadata
KIR <- AddMetaData(KIR, demux_data$Barcodes_classification, col.name = 'barcode')

# Add barcode matrix as assay
# Create Barcode assay ?
KIR[["Barcodes"]] <-  CreateAssayObject(counts = BC_raw)

# same as above - table from demux_data which is now in KIR object
table(KIR$barcode)

# Remove barcode and hashing rows from ADT data
phenotype.list <- rownames(KIR[["ADT"]])[-26:0]

# Extract barcode data from Antibody Capture matrix
phenotype <- as.matrix(GetAssayData(object = KIR[["ADT"]], slot = "counts"))[phenotype.list,]

# Create assay object while overriding old ADT assay
KIR[["ADT"]] <- CreateAssayObject(counts = phenotype)

# Create an object of KIR without doublets:
KIR.subset <- KIR %>% subset(HTO_classification.global == "Singlet") #Now we have 4655 instead of 6360 GEMs

# Add conditions for each sample as metadata:
KIR.subset$condition <- recode(KIR.subset$HTO_classification,
                               "S8"="KIR+_MM+",  "S2"="MM+",
                               "S4"="KIR+_MM+", "S5"="control",
                               "S6"="TCR_specific", "S3"="MM+",
                               "S10"="control", "S1"="TCR_specific",
                               "S7"="TCR_specific", "S9"="KIR+_MM+")
KIR$condition <- recode(KIR$HTO_classification,
                        "S8"="KIR+_MM+",  "S2"="MM+",
                        "S4"="KIR+_MM+", "S5"="control",
                        "S6"="TCR_specific", "S3"="MM+",
                        "S10"="control", "S1"="TCR_specific",
                        "S7"="TCR_specific", "S9"="KIR+_MM+")

# Venn diagrams for data overview
pMHCdata <-  head(KIRdata[["Antibody Capture"]],16)
i <- (colSums(pMHCdata, na.rm=T) > 10) #filter out columns with less than 10
pMHCdata <- pMHCdata[, i]

ADTdata <- tail(KIRdata[["Antibody Capture"]],137)
p <- (colSums(ADTdata, na.rm=T) > 300) #filter out colSums under 300
ADTdata <- ADTdata[, p]

HTOdata <- KIR@assays[["HTO"]]@counts

RNAdata <- KIRdata[["Gene Expression"]]

list_venn <- list(RNA = colnames(RNAdata),
                  ADT = colnames(ADTdata),
                  HTO = colnames(HTOdata),
                  pMHC = colnames(pMHCdata),
                  TCR = vdj_data$barcode)
#percentage
ggVennDiagram(list_venn, label_alpha = 0, label = "percent") +
  scale_fill_gradient(low="lightskyblue",high = "yellow")

#counts
venn(list_venn, ilab=TRUE, zcolor = "style")



############################### pMHC Barcode filtering ####################################
exp_pairs <- c('BC2_BC3', 'BC4_BC5')
KIR.subset$bcexp_pair_or_sin <- ifelse(!grepl('_', KIR.subset$barcode) | KIR.subset$barcode %in% exp_pairs,
 T, F)
KIR.subset$bcexp_pair_or_sin[is.na(KIR.subset$barcode)] <- NA

KIR.subset$barcode[(KIR.subset$barcode == "BC4_BC6") & (KIR.subset@assays$Barcodes@counts["BC4"<2])]


KIR.subset$barcode[(KIR.subset$HTO_classification == "S9") & 
                     (KIR.subset$barcode=="BC14_BC2") &
                     (!is.na(KIR.subset$barcode))] <-"BC14"

KIR.subset$barcode[(KIR.subset$HTO_classification == "S9") & 
                     (KIR.subset$barcode=="BC14_BC7") &
                     (!is.na(KIR.subset$barcode))] <-"BC14"

KIR.subset$barcode[(KIR.subset$HTO_classification == "S9") & 
                     (KIR.subset$barcode=="BC14_BC6") &
                     (!is.na(KIR.subset$barcode))] <-"BC14"

KIR.subset$barcode[(KIR.subset$HTO_classification == "S9") & 
                     (KIR.subset$barcode=="BC14_BC5") &
                     (!is.na(KIR.subset$barcode))] <-"BC14"


# Make barcode count plot
KIR.subset@meta.data %>% 
  filter(bcexp_pair_or_sin) %>%
  dplyr::select(c(hash.ID, barcode, buffy_coat)) %>% 
  filter(barcode > 10) %>%  # remove those with less than 10 barcodes
  filter(!is.na(hash.ID)) %>% 
  mutate(buffy_coat = factor(buffy_coat)) %>% 
  mutate(barcode = factor(barcode)) %>% 
  group_by(buffy_coat, barcode, hash.ID) %>% 
  summarise(barcodecount=n()) %>% 
  ggplot(aes(x=hash.ID, y=barcodecount, fill=barcode)) +
  xlab("Hashing") + ylab("Barcode count")+
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  geom_text(aes(label = barcode), position = position_stack(vjust = 0.5), size=3 ) +
  theme_minimal()+
  facet_grid(~buffy_coat, scales = "free")



#################### Standard pre-processing workflow ##################

# 1)  QC and Selection of Cells for Analysis ----------------------------------------------------

# Save percentage of mitochondrial genes to metadata:
KIR.subset[["percent.mt"]] <- PercentageFeatureSet(KIR.subset, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(KIR.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter plots of QC metrics
plot1 <- FeatureScatter(KIR.subset, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  xlab("RNA count") + ylab("Percentage mitochondrial genes")
                              
plot2 <- FeatureScatter(KIR.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_smooth(method = "lm") 
  xlab("RNA count") + ylab("Gene count")
                              
plot1 + plot2



# 2)  Filtering out low quality cells ------------------------------------------------------------
# we filter out cells that have unique feature counts over 2,500 or less than 200
# we filter out cells that have more than 5% mitochondrial counts
KIR.subset <-  subset(KIR.subset, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 
KIR.subset #4143 instead of 4655 


# 3)  Normalizing the data -------------------------------------------------------------------------
# lognormalized and multiplied by 10000, these are default
KIR.subset <- NormalizeData(KIR.subset, normalization.method = "LogNormalize", scale.factor = 10000)


# Remove non-coding genes from KIR.subset object
genes <- c(KIR@assays[["RNA"]]@data@Dimnames[[1]])
require(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "https://useast.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
all_coding_genes <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol",
    "gene_biotype"),
  filter = "hgnc_symbol",
  values = genes,
  uniqueRows=TRUE,
)
all_coding_genes <- all_coding_genes  %>% subset(gene_biotype=="protein_coding")

# Create a seurat object with the clean data:
KIR.subset@assays[["RNA"]]@counts <- KIR.subset@assays[["RNA"]]@counts[all_coding_genes$hgnc_symbol ,]

# Check how many genes are left in KIR.subset (should be around 20,000 protein-coding genes in humans)
nrow(KIR.subset@assays[["RNA"]]@counts)
# 18216 genes left 



# 4)  Identification of highly variable features -----------------------------------------------------
# We wish to find genes that have high cell to cell variation
KIR.subset <- FindVariableFeatures(KIR.subset, selection.method = "vst", nfeatures = 2000) 

# We then remove TCR genes found by the regex:
KIR.subset@assays[["RNA"]]@var.features <- KIR.subset@assays[["RNA"]]@var.features[!grepl(paste0(
                                            '^IG[H,K,L][A-K, M-Z].*|^TR[A,B][D,J,V].*'
                                            , collapse = "|"), KIR.subset@assays[["RNA"]]@var.features)]

# Identify the 10 most highly variable genes in both objects
top10.subset <- head(VariableFeatures(KIR.subset), 10)
top10.subset

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(KIR.subset)
LabelPoints(plot = plot1, points = top10, repel = TRUE) #the red points are the features found by FindVariableFeatures



# 5)  Scaling the Data ----------------------------------------------------------------------------
# we need to scale the data as there are many sources of unwanted variation
# this could be technical noice (e.g. batch effects) or due to biological sources (e.g. differences in cell cycles)
# we want to account for these variations to prevent cells from clustering together due to a variation 
all.genes.subset <- rownames(KIR.subset)
KIR.subset <- ScaleData(KIR.subset, features = all.genes.subset)



# 6)  Perform Linear Dimensional Reduction (PCA) ---------------------------------------------------
#this is done in order to identify the sources of heterogenity in the data set
#by default it only considers the 2000 variable features 
KIR.subset <- RunPCA(KIR.subset, features = VariableFeatures(object = KIR.subset))

# Examine and visualize PCA results a few different ways
DimPlot(KIR.subset, reduction = "pca")

# Heatmaps can be used to identify which PCs we should use in downstream analysis
#cells and features are ordered according to PCA scores
#only 500 cells are included in both ends of the spectrum to increase speed
DimHeatmap(KIR.subset, dims = 1:15, cells = 500, balanced = TRUE) #15 PCs

# Choosing only significant PCs which capture the majority of the signal in downstream analysis
ElbowPlot(KIR.subset)



# 7)  Cluster the Cells -----------------------------------------------------------------------
#compute nearest neighbours (KNN) and SNN
KIR.subset <- FindNeighbors(KIR.subset, dims = 1:15)

# Setting the resolution
KIR.subset <- FindClusters(KIR.subset, resolution = c(0.1, 0.2,0.3, 0.4, 0.5, 0.7, 1))

# Vizualise how many clusters there are for each resolution
DimPlot(KIR.subset, group.by="RNA_snn_res.0.5", label = TRUE)
Idents(KIR.subset) <- "RNA_snn_res.0.5"

# Create UMAPs
KIR.subset <- RunUMAP(KIR.subset, dims = 1:15)
DimPlot(KIR.subset,reduction = "umap", label=T, label.size = 6)

# UMAP by condition:
DimPlot(KIR.subset, reduction = "umap", label = T, group.by = "RNA_snn_res.0.5", split.by = "condition") +
  ggtitle("UMAP by condition")

# UMAP by sample
DimPlot(KIR.subset, group.by = "hash.ID", label = F) + 
  ggtitle('All samples', subtitle = "coloured based on series") + 
  facet_grid("hash.ID")
  theme(plot.title = element_text(hjust = 0,size = 15),
        plot.subtitle = element_text(hjust = 0, size = 15))

# UMAP by subject
DimPlot(KIR.subset, group.by = "buffy_coat", label = F) + 
  ggtitle('UMAP grouped by subject') + 
  theme(plot.title = element_text(hjust = 0,size = 15),
        plot.subtitle = element_text(hjust = 0, size = 15))

# UMAP by barcode
DimPlot(KIR.subset, group.by = "barcode", label = F) + 
  ggtitle('UMAP grouped by barcode') + 
  theme(plot.title = element_text(hjust = 0,size = 15),
        plot.subtitle = element_text(hjust = 0, size = 15))



# 8) Differential expression analysis (find biomarkers) -----------------------------------------
DefaultAssay(KIR.subset) <- "RNA"
Idents(KIR.subset) <- KIR.subset$condition

# Find markers for every condition compared to all remaining cells, report only the positive ones
KIR.markers.subset <- FindAllMarkers(KIR.subset, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

# Find markerrs for conditions compared to other specific conditions
# 1) KIR+_MM+ vs TCR specific
# 2) MM+ vs KIR+_MM+
# 3) TCR specific vs KIR+_MM+, MM+, control
KIR.markers1.subset <- FindMarkers(KIR.subset, 
                                   ident.1 = "KIR+_MM+", ident.2 = "TCR_specific", 
                                   only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  mutate(p_val_adj = p.adjust(p_val, method = 'BH'))
KIR.markers2.subset <- FindMarkers(KIR.subset, ident.1 = "MM+", ident.2 = "KIR+_MM+", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)  %>%
  mutate(p_val_adj = p.adjust(p_val, method = 'BH'))
KIR.markers3.subset <- FindMarkers(KIR.subset, ident.1 = "TCR_specific", ident.2 = c("KIR+_MM+","MM+","control") ,only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  mutate(p_val_adj = p.adjust(p_val, method = 'BH'))


# Remove TCR genes from the markers
KIR.markers1.subset <- KIR.markers1.subset %>% 
  filter(!str_detect(row.names(KIR.markers1.subset), '^IG[H,K,L][A-K, M-Z].*|^TR[A,B,G,D][D,J,V,C].*'))
KIR.markers2.subset <- KIR.markers2.subset %>% 
  filter(!str_detect(row.names(KIR.markers2.subset), '^IG[H,K,L][A-K, M-Z].*|^TR[A,B,G,D][D,J,V,C].*'))
KIR.markers3.subset <- KIR.markers3.subset %>% 
  filter(!str_detect(row.names(KIR.markers3.subset), '^IG[H,K,L][A-K, M-Z].*|^TR[A,B,G,D][D,J,V,C].*'))

# Remove ribosomal gene
KIR.markers1.subset <- KIR.markers1.subset %>% 
  filter(!str_detect(row.names(KIR.markers1.subset), "RP11-277P12.6"))

# Gene expression shown for each condition for top 10 significant genes orderd by log2FC:
feature.list1 <- KIR.markers1.subset %>% 
  filter(p_val_adj < .05) %>%
  slice_max(n = 5, order_by=avg_log2FC) %>%
  row.names()

feature.list2 <- KIR.markers2.subset %>% 
  filter(p_val_adj < .05) %>%
  slice_max(n = 5, order_by=avg_log2FC) %>%
  row.names()

feature.list3 <- KIR.markers3.subset %>% 
  filter(p_val_adj < .05) %>%
  slice_max(n = 5, order_by=avg_log2FC) %>%
  row.names()

# Create heatmap showing gene expressions 
KIR.subset$condition <- factor(KIR.subset$condition, levels = c('KIR+_MM+','control', 'MM+', 'TCR_specific'))
DoHeatmap(KIR.subset, group.by = 'condition', slot = 'scale.data',
          features = c(feature.list1, "KIR3DL1", "B2M","HLA-C", feature.list2, feature.list3), angle=0, size=4, hjust = 0.5)+ NoLegend() + 
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(text = element_text(size = 15, vjust = 0.5)) +   
  guides(color='none') 

#Show the significant markers
KIR.markers.subset %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) 

# Visualizes feature expression on a tSNE or PCA plot
DefaultAssay(KIR.subset) = "RNA"
FeaturePlot(KIR.subset, features = c("KIR3DL1", "KIR2DL1", "XCL2", "GNLY"))

Idents(KIR.subset) <- KIR.subset$buffy_coat
FeaturePlot(KIR.subset, features = c("KIR3DL1", "KIR2DL1"), split.by = "condition", label = T)

# GD T cells test
FeaturePlot(KIR.subset, c('TRDC', "TRDV1")) 



##################### Differential expression analysis of surface markers ######################
# Normalize and scale the data
KIR.subset <- KIR.subset %>% 
  NormalizeData(normalization.method='CLR', assay ='ADT') %>% 
  ScaleData(assay='ADT')

# Find markers of all and chosen comparisons
# 1) KIR+_MM+ vs TCR specific
# 2) MM+ vs KIR+_MM+
# 3) TCR specific vs KIR+_MM+, MM+, control
Idents(KIR.subset) <- KIR.subset$condition 
ADT.markers <- FindAllMarkers(KIR.subset,assay = "ADT" ,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

ADT.markers1 <- FindMarkers(KIR.subset, assay = "ADT", ident.1 = "KIR+_MM+", ident.2 = "TCR_specific" ,only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.1) %>%
  mutate(p_val_adj = p.adjust(p_val, method = 'BH'))
ADT.markers2 <- FindMarkers(KIR.subset, assay = "ADT", ident.1 = "MM+", ident.2 = "KIR+_MM+" ,only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.1) %>%
  mutate(p_val_adj = p.adjust(p_val, method = 'BH'))
ADT.markers3 <- FindMarkers(KIR.subset, assay = "ADT",ident.1 = "TCR_specific", ident.2 = c("KIR+_MM+","MM+","control") ,only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.1) %>%
  mutate(p_val_adj = p.adjust(p_val, method = 'BH'))

# Remove TCR genes
ADT.markers3 <- ADT.markers3 %>% 
  filter(!str_detect(row.names(ADT.markers3), '^IG[H,K,L][A-K, M-Z].*|^TR[A,B,G,D][D,J,V,C].*'))

# Gene expression shown for each condition for top 10 significant genes orderd by log2FC:
ADT.feature.list1 <- ADT.markers1 %>% 
  filter(p_val_adj < .05) %>%
  slice_max(n = 5, order_by=avg_log2FC) %>%
  row.names()

ADT.feature.list2 <- ADT.markers2 %>% 
  filter(p_val_adj < .05) %>%
  slice_max(n = 5, order_by=avg_log2FC) %>%
  row.names()

ADT.feature.list3 <- ADT.markers3 %>% 
  filter(p_val_adj < .05) %>%
  slice_max(n = 5, order_by=avg_log2FC) %>%
  row.names()

# Create heatmap showing surface antigen expressions 
KIR.subset$condition <- factor(KIR.subset$condition, levels = c('KIR+_MM+','control','MM+','TCR_specific' ))
DoHeatmap(KIR.subset, assay = "ADT", group.by = 'condition',
          features = c(ADT.feature.list1, ADT.feature.list2,ADT.feature.list3), angle=0, size=4, hjust = 0.5)+ NoLegend() + 
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(text = element_text(size = 15, vjust = 0.5)) +   
  guides(color='none')



############################ Clonotype Analysis #############################3
# Clonotype plot with barcodes as fill
KIR.subset$hash.ID <- factor(KIR.subset$hash.ID, levels=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10"))
Idents(KIR.subset) <- "hash.ID"

ID <- KIR.subset$HTO_classification
clonotype <- KIR.subset$TRB_clone
barcode <- KIR.subset$barcode

f <- as.data.frame(cbind(ID,clonotype,barcode))
f <- f[!is.na(f$clonotype),]
f <- arrange(f, ID)

clono_levels <- f %>% pull(clonotype) %>% unique()

f <- mutate(f, clonotype = factor(clonotype,levels = clono_levels))

table(KIR.subset$HTO_classification.global)

ggplot(f, aes(x = clonotype, y = ID)) +
  geom_count(aes(color = barcode)) + 
  theme(axis.text.x = element_blank()) + xlab("clonotype") + ylab("sample") +
  theme_minimal() +
  scale_y_discrete(limits=rev) +
  facet_grid(buffy_coat)



# Investigation of top 10 clonotypes -----------------------------------------------------------

# Top 10 clones for stacked barplots:
#Buffy coat 357
top10_b_BC357 <- KIR.subset@meta.data %>% subset(buffy_coat == "BC357")
top10_b_BC357 <- sort(table(top10_b_BC357$TRB_clone), decreasing = T)[1:10]
topb_clones_BC357 <- paste0("clone=", names(top10_b_BC357), ", size=", top10_b_BC357)
names(topb_clones_BC357) <- names(top10_b_BC357)
topb_clones_BC357 <- as.data.frame(topb_clones_BC357) %>% 
  rownames_to_column("TRB_clone") 

all(Cells(KIR.subset) == rownames(KIR.subset@meta.data))

KIR.subset@meta.data <- KIR.subset@meta.data %>%
  mutate(TRB_clone = as.character(TRB_clone)) %>%
  left_join(topb_clones_BC357, by= "TRB_clone")
rownames(KIR.subset@meta.data) <- Cells(KIR.subset) 

#Buffy coat 351
top10_b_BC351 <- KIR.subset@meta.data %>% subset(buffy_coat == "BC351")
top10_b_BC351 <- sort(table(top10_b_BC351$TRB_clone), decreasing = T)[1:10]
topb_clones_BC351 <- paste0("clone=", names(top10_b_BC351), ", size=", top10_b_BC351)
names(topb_clones_BC351) <- names(top10_b_BC351)
topb_clones_BC351 <- as.data.frame(topb_clones_BC351) %>% 
  rownames_to_column("TRB_clone") 

all(Cells(KIR.subset) == rownames(KIR.subset@meta.data))

KIR.subset@meta.data <- KIR.subset@meta.data %>%
  mutate(TRB_clone = as.character(TRB_clone)) %>%
  left_join(topb_clones_BC351, by= "TRB_clone")
rownames(KIR.subset@meta.data) <- Cells(KIR.subset) 


# UMAP showing the clones:
DimPlot(KIR.subset, label = T, label.size = 3, group.by = "topb_clones", cols=c("#00B0F6" ,"#FF62BC", "#00BF7D", "#00BFC4", "#39B600", "#D89000", "#F8766D","#9590FF", "#A3A500" ,"#E76BF3"), na.value = "grey90")

# Clonotype plot 
IDs <- data.frame(hash.ID=KIR.subset$HTO_classification) %>% tibble::rownames_to_column("barcode") %>% as_tibble()
buffycoats <- data.frame(buffycoat =KIR.subset$buffy_coat ) %>% tibble::rownames_to_column("barcode") %>% as_tibble()
conditions <- data.frame(condition =KIR.subset$condition ) %>% tibble::rownames_to_column("barcode") %>% as_tibble()
clonotypes <- vdj_data %>% select(barcode,TRB_clone) %>% as_tibble()

df_clonotype <- left_join(clonotypes, IDs,by='barcode') %>% 
  drop_na(hash.ID) %>% 
  arrange(hash.ID) %>% 
  filter(hash.ID %in% c('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10'))

df_clonotype_BC351 <- full_join(df_clonotype, buffycoats ,by='barcode') %>% 
  drop_na(buffycoat) %>% 
  arrange(buffycoat) %>% 
  filter(buffycoat %in% c('BC351'))

df_clonotype_BC357 <- full_join(df_clonotype, buffycoats ,by='barcode') %>% 
  drop_na(buffycoat) %>% 
  arrange(buffycoat) %>% 
  filter(buffycoat %in% c('BC357'))

df_clonotype_BC351 <- full_join(df_clonotype, conditions ,by='barcode') %>% 
  drop_na(condition) %>% 
  arrange(condition) 

df_clonotype_BC357 <- full_join(df_clonotype, conditions ,by='barcode') %>% 
  drop_na(condition) %>% 
  arrange(condition) 

clono_levels_BC351 <- df_clonotype_BC351 %>%
  pull(TRB_clone) %>%
  unique()

clono_levels_BC357 <- df_clonotype_BC357 %>%
  pull(TRB_clone) %>%
  unique()

df_clonotype_BC351 <- mutate(df_clonotype_BC351, clonotype = factor(TRB_clone, levels = clono_levels_BC351))
df_clonotype_BC357 <- mutate(df_clonotype_BC357, clonotype = factor(TRB_clone, levels = clono_levels_BC357))

df_clonotype_BC351 <- filter(df_clonotype_BC351, TRB_clone %in% names(top10_b_BC351) )
df_clonotype_BC357 <- filter(df_clonotype_BC357, TRB_clone %in% names(top10_b_BC357) )

df_clonotype_BC357$hash.ID <- factor(df_clonotype_BC357$hash.ID, 
                      levels = c("S1","S2","S3","S4","S5",
                                 "S6","S7","S8","S9","S10"))

# Create plots with point size and colour relfecting counts 
ggplot(df_clonotype_BC351, aes(x = clonotype, y = hash.ID)) +
  geom_count(aes(color = ..n.., size = ..n..), alpha=0.8) + 
  scale_y_discrete(limits=rev) +
  scale_size_area(max_size = 10, breaks=c(1,2,3,4,5,6,7)) + 
  scale_color_viridis(option = "D", begin = 0, end = .75, 
                      breaks = c(11,2,3,4,5,6,7), labels = c(1,2,3,4,5,6,7)) + 
  theme_minimal() + 
  theme(text = element_text(size=20), strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15)) + 
  labs(x="beta clone ID",y="",size='GEMs',color='GEMs') +
  ggtitle("Buffycoat 351")

ggplot(df_clonotype_BC357, aes(x = clonotype, y = hash.ID)) +
  geom_count(aes(color = ..n.., size = ..n..), alpha=0.8) + 
  scale_y_discrete(limits=rev) +
  scale_size_area(max_size = 10, breaks=c(1,10,50,100,150,200,250,300,350,400)) + 
  scale_color_viridis(option = "D", begin = 0, end = .75, 
                      breaks = c(1,10,50,100,150,200,250,300,350,400), labels = c(1,10,50,100,150,200,250,300,350,400)) + 
  theme_minimal() + 
  theme(text = element_text(size=25), strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15)) + 
  labs(x="beta clone ID",y="",size='GEMs',color='GEMs') +
  ggtitle("Buffycoat 357")


# Create a similar plot with percentage as size instead:
df_clonotype_percent <- df_clonotype_BC351%>%
  group_by(hash.ID, clonotype) %>% mutate(GEM_count=n()) %>%
  ungroup() %>% mutate(GEM_sum=sum(GEM_count)) %>%
  mutate(percent=(GEM_count/GEM_sum)*100) #%>%

ggplot(df_clonotype_percent,
             aes(x = clonotype, y = hash.ID, size=percent, color=percent)) +
  geom_point(alpha=0.8) +
  scale_y_discrete(limits=rev) +
  scale_size_area(max_size = 25, breaks=c(0.1,1, 5, 10,20,30)) +
  scale_color_viridis(option = "D", begin = 0, end = .75, direction = -1,
                      breaks = c(10,20,30,40), labels = c(10,20,30,40)) +
  theme_bw() +
  theme(text = element_text(size=20), strip.background = element_blank()) +
  theme(legend.key.size = unit(.5, 'cm'), #change legend key size
        legend.key.height = unit(.5, 'cm'), #change legend key height
        legend.key.width = unit(.5, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) +#change legend text font size
  labs(x="clonotype ID",y="",size='percent',color='percent') +
  ggtitle("Buffy coat 351") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# Sankey plots -----------------------------------------------------------
# NB - Made in python 

#Alpha vs Beta-chain - use VDJ_table.csv

# Beta-chain clones, conditions and barcodes table to load in python
top20_b <- sort(table(KIR.subset$TRB_clone), decreasing = T)[1:20]
top20_beta_clones <- paste0(names(top20_b))
names(top20_beta_clones) <- names(top20_b)
top20_beta_clones <- as.data.frame(top20_beta_clones) %>% 
  rownames_to_column("TRB_clone") 

all(Cells(KIR.subset) == rownames(KIR.subset@meta.data))

KIR.subset@meta.data <- KIR.subset@meta.data %>%
  mutate(TRB_clone = as.character(TRB_clone)) %>%
  left_join(top20_beta_clones, by= "TRB_clone")
rownames(KIR.subset@meta.data) <- Cells(KIR.subset) 

sankey_df <- data.frame (HTO_maxID  = KIR.subset$HTO_maxID  ,
                         barcodes = KIR.subset$barcode  ,
                         condition = KIR.subset$condition  ,
                         clonotype =  KIR.subset$top20_beta_clones  
)

# Save dataframes for import in Python:
write.csv(sankey_df,"/Users/simone/Desktop/MASTER/sankey_data.csv", row.names = T)
write.csv(sankey_df,"/Volumes/T-cells-and-cancer/SRH group/Group members/Simone/data/metadata/sankey_data.csv", row.names = T)



############################ Doublets coming from hashing-antibodies #############################
# Find the most frequent hashing-antibody in thee doublets 
KIR$HTO_classification[KIR$HTO_classification.global == "Doublet"]

a <- separate(as_tibble(KIR$HTO_classification[KIR$HTO_classification.global == "Doublet"]), value , c("N1","N2"), sep="_")
b <- a %>% count(N1)  
c <- a %>% count(N2)  
d <- full_join(b,c, by = c("N1"="N2")) %>% replace(is.na(.), 0) 

d <- d %>% as_tibble() %>% 
  mutate(sum = rowSums(across(where(is.numeric))))
# add frequecy
k <- (d$sum/2966)*100
d$freq <- k

# Make barplot with count and frequency
ggplot(data=d, aes(x=N1, y=sum)) +
  geom_bar(stat="identity", fill=c("forestgreen"))+
  xlab("Antibodies") + ylab("Count")+
  theme_minimal()

ggplot(data=d, aes(x=N1, y=freq)) +
  geom_bar(stat="identity", fill=c("forestgreen"))+
  xlab("Antibodies") + ylab("Frequency in percent")+
  theme_minimal()



############################ Gene Enrichment Analysis #############################
Idents(KIR.subset) <- KIR.subset$condition

# Keep only the significant genes and divide into conditions
# Use already found marker genes within each comparison:
# 1) KIR+_MM+ vs TCR specific
# 2) MM+ vs KIR+_MM+
# 3) TCR specific vs KIR+_MM+, MM+, control
markers_sig1 <- subset(KIR.markers1.subset, p_val_adj < 0.05)
markers_sig2 <- subset(KIR.markers2.subset, p_val_adj < 0.05)
markers_sig3 <- subset(KIR.markers3.subset, p_val_adj < 0.05)

markers1_KIR <- markers_sig1 %>% filter(str_detect(cluster, "KIR+_MM+"))
markers1_TCR_specific <- markers_sig1 %>% filter(str_detect(cluster, "TCR_specific"))
markers1_control <- markers_sig1 %>% filter(str_detect(cluster, "control"))
markers1_MM <- markers_sig1 %>% filter(str_detect(cluster, "MM+"))

markers2_KIR <- markers_sig2 %>% filter(str_detect(cluster, "KIR+_MM+"))
markers2_TCR_specific <- markers_sig2 %>% filter(str_detect(cluster, "TCR_specific"))
markers2_control <- markers_sig2 %>% filter(str_detect(cluster, "control"))
markers2_MM <- markers_sig2 %>% filter(str_detect(cluster, "MM+"))

markers3_KIR<- markers_sig3 %>% filter(str_detect(cluster, "KIR+_MM+"))
markers3_TCR_specific <- markers_sig3 %>% filter(str_detect(cluster, "TCR_specific"))
markers3_control <- markers_sig3 %>% filter(str_detect(cluster, "control"))
markers3_MM<- markers_sig3 %>% filter(str_detect(cluster, "MM+"))

# Get the significant up-regulated genes
up <- subset(markers_sig, avg_log2FC > 0)
up_KIR <- subset(markers_KIR, avg_log2FC > 0)
up_TCR_specific <- subset(markers_TCR_specific, avg_log2FC > 0)
up_control <- subset(markers_control, avg_log2FC > 0)
up_MM <- subset(markers_MM, avg_log2FC > 0)


# Multi analysis
multi_gp_KIR <- gost(list("up-regulated - KIR and multimer positive" = row.names(up_KIR),
                     sources = c("GO:MF", "GO:BP", "GO:CC", "KEGG","REAC", "HPA")))
plot_KIR <- gostplot(multi_gp_KIR, interactive = T)
plot_KIR


multi_gp_TCR_specific <- gost(list("up-regulated - TCR_specific" = row.names(up_TCR_specific),
                      sources = c("GO:MF", "GO:BP", "GO:CC", "KEGG","REAC", "HPA")))
plot_TCR_specific<- gostplot(multi_gp_TCR_specific, interactive = T)
plot_TCR_specific


multi_gp_control <- gost(list("up-regulated - control" = row.names(up_control),
                      sources = c("GO:MF", "GO:BP", "GO:CC", "KEGG","REAC", "HPA")))
plot_control<- gostplot(multi_gp_control, interactive = T)
plot_control


multi_gp_MM <- gost(list("up-regulated - Multimer positive" = row.names(up_MM),
                      sources = c("GO:MF", "GO:BP", "GO:CC", "KEGG","REAC", "HPA")))
plot_MM<- gostplot(multi_gp_MM, interactive = T)
plot_MM


################################ Barcode threshold checks ########################################
# Scatterplot with counts for each barcode 
# Create a random colour palette to distinguish all the barcodes
n <- 26
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

# Use the barcode assat
DefaultAssay(KIR.subset) <- "Barcodes"

# Make scatterplot
FeatureScatter(KIR.subset, feature1 = "BC1", feature2 = "BC6", group.by = "barcode") +
  scale_color_manual(values=sample(col_vector, n))


mat <- KIR.subset@assays$Barcodes@counts[,KIR.subset$HTO_classification == "S7"]
barcodes <- KIR.subset@meta.data$barcode[KIR.subset$HTO_classification == "S7"]

qplot(mat['BC6',], mat['BC7',], color=barcodes) +
  scale_color_manual(values=sample(col_vector, n)) +
  xlab("BC6") + ylab("BC7")+
  theme_bw() 



#-----------------------------------------------------------------------------------------------------
save.image('/Users/simone/Desktop/KIR_Simone.RData')
