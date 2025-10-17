#reference 
#https://nbisweden.github.io/workshop-scRNAseq/labs/seurat/seurat_01_qc.html#meta-qc_doublet
#https://support.parsebiosciences.com/hc/en-us/articles/360053078092-Seurat-Tutorial-65k-PBMCs



library(Seurat)
library(dplyr)
library(tibble)
library(Matrix)
library(ggplot2)
library(DoubletFinder)
library(devtools)

rm(list = ls())


data_path <- "C:/workspace/Seurat/mclover3/data/"
fig_path <- "C:/workspace/Seurat/mclover3/figure/"

# Convenience functions
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}

DGE_folder <- "C:/workspace/Seurat/mclover3/mclover3_D10_12_2024-07-09/DGE_filtered/"
#DGE_folder <- "C:/workspace/Seurat/mclover3/mclover3_target_well_2024-06-14/DGE_filtered/"
#DGE_folder <- "C:/Users/jhuang45/Downloads/parse_splitseq_pipleline/all-sample/DGE_filtered/"
#DGE_folder <- "C:/Users/jhuang45/Downloads/DGE_filtered_20231213/"
# split-pipe versions older than 1.1.0 used "DGE.mtx"
mat <- readMM(paste0(DGE_folder, "count_matrix.mtx"))

cell_meta <- read.delim(paste0(DGE_folder, "cell_metadata.csv"),
                        stringsAsFactor = FALSE, sep = ",")
genes <- read.delim(paste0(DGE_folder, "all_genes.csv"),
                    stringsAsFactor = FALSE, sep = ",")

cell_meta$bc_wells <- make.unique(cell_meta$bc_wells, sep = "_dup")
rownames(cell_meta) <- cell_meta$bc_wells
genes$gene_name <- make.unique(genes$gene_name, sep = "_dup")

# Setting column and rownames to expression matrix
colnames(mat) <- genes$gene_name
rownames(mat) <- rownames(cell_meta)
mat_t <- t(mat)

# Remove empty rownames, if they exist
mat_t <- mat_t[(rownames(mat_t) != ""),]

# Seurat version 5 or greater uses "min.features" instead of "min.genes"
all_sample <- CreateSeuratObject(mat_t, min.features = 300, min.cells = 5, meta.data = cell_meta)

# Setting our initial cell class to a single type, this will changer after clustering. 
all_sample@meta.data$orig.ident <- factor(rep("all_sample", nrow(all_sample@meta.data)))
Idents(all_sample) <- all_sample@meta.data$orig.ident

#SaveObject(all_sample, "seurat_obj_before_QC")
#all_sample <- ReadObject("seurat_obj_before_QC")


#Cell quality control
#all_sample[["percent.mt"]] <- PercentageFeatureSet(all_sample, pattern = "^MT-")
#all_sample[["percent.mt"]] <- PercentageFeatureSet(all_sample, pattern = "^Mt-|^AY172581")
all_sample[["percent.mt"]] <- PercentageFeatureSet(all_sample, pattern = "^(ATP6|ATP8|COX[1-3]|CYTB|mt-Rnr[12]|ND[1-6L]?|Trn[a-y](1|2)?)$")
head(all_sample)
#alldata <- PercentageFeatureSet(all_sample, "^MT-", col.name = "percent_mito")
#head(alldata, 10)


plot <- VlnPlot(all_sample, pt.size = 0.10,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SaveFigure(plot, "vln_QC", width = 12, height = 6)

plot1 <- FeatureScatter(all_sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all_sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
SaveFigure((plot1 + plot2),"scatter_QC", width = 12, height = 6, res = 200)
plot
plot1
plot2
mean(all_sample$percent.mt)
sd(all_sample$percent.mt)
# Perform the filtering, change nFeature_RNA & nCount_RNA  & percent.mt according plot
all_sample <- subset(all_sample, subset = nFeature_RNA < 2000 & nCount_RNA < 3000 & percent.mt < 3)
dim(all_sample)



#Normalizing the data
all_sample <- NormalizeData(all_sample, normalization.method = "LogNormalize", scale.factor = 10000)

#all_sample <- CellCycleScoring(object = all_sample, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)

#Identification of highly variable features
all_sample <- FindVariableFeatures(all_sample, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(all_sample), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(all_sample)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
SaveFigure((plot1 + plot2), "var_features", width = 12, height = 6)
plot1+plot2
plot1
plot2
#Scaling the data
all_sample <- ScaleData(all_sample)
#Perform linear dimensional reduction
all_sample <- RunPCA(all_sample )

#SaveObject(all_sample, "seurat_obj_after_PCA")

#all_sample <- ReadObject("seurat_obj_after_PCA")
head(all_sample)
# Examine and visualize PCA results a few different ways
print(all_sample[["pca"]], dims = 1:5, nfeatures = 5)
plot <- VizDimLoadings(all_sample, dims = 1:2, reduction = "pca")
SaveFigure(plot, "viz_PCA_loadings", width = 10, height = 8)
plot





SaveObject(all_sample, "seurat_obj_before_cluster_mclover3_target_well")
all_sample <- ReadObject("seurat_obj_before_cluster_mclover3_target_well")
#Cluster the cells
all_sample <- FindNeighbors(all_sample, dims = 1:30)
all_sample <- FindClusters(all_sample, resolution = 1.5) #change resolution to get different clusters
all_sample <- BuildClusterTree(all_sample, reorder = TRUE, reorder.numeric = TRUE)

#SaveObject(all_sample, "seurat_obj_after_cluster_illumina")
#all_sample <- ReadObject("seurat_obj_after_cluster_illumina")

cluster <- all_sample@meta.data %>% 
  rownames_to_column("barcodes") %>%
  select(barcodes, tree.ident)
colnames(cluster) <- c('barcode', 'cluster.mclover3')
write.csv(cluster, "cluster_mclover3.csv",row.names = FALSE)

cluster_well <- all_sample@meta.data %>% 
  rownames_to_column("barcodes") %>%
  select(barcodes, tree.ident,bc_wells,bc1_well,bc2_well,bc3_well)
colnames(cluster_well) <- c('barcode', 'cluster.mclover3')
write.csv(cluster_well, "cluster_mclover3.csv",row.names = FALSE)

a<-cluster_well[c("barcodes","bc1_well")]
unique(a)

all_sample <- RunUMAP(all_sample, dims = 1:30)
pdf(file="umap_mclover3_resolution_50.pdf", width=4.4,height=4.76)
DimPlot(all_sample, reduction = "umap", label = TRUE)+ NoLegend()
VlnPlot(all_sample, features = c("Drd1", "Drd2","Pde4b"), group.by = "tree.ident")
dev.off()
DimPlot(all_sample, reduction = "umap", label = TRUE) + NoAxes() + NoLegend()


#Differential gene expression (finding cluster markers)
all_markers <- FindAllMarkers(all_sample, min.pct = 0.25, logfc.threshold = 0.25)
all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
VlnPlot(all_sample, features = c("Drd1", "Drd2","Pde4b"), group.by = "tree.ident")
VlnPlot(all_sample, features = c("Drd1", "Drd2"), slot = "counts", log = TRUE)

#Visualizing the top n genes per cluster
top1 <- all_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top5 <- all_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top20 <- all_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top30 <- all_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

to_plot <- unique(top5$gene)
to_plot
DotPlot(all_sample, features = to_plot, group.by = "tree.ident") + coord_flip()
to_plot_top20<-unique(top20$gene)
to_plot_top20
to_plot_top30<-unique(top30$gene)
to_plot_top30
#to_plot_top50<-unique(top50$gene)

library("stringr")

gabaergic_interneuron<-str_to_title(c("ADAM22", "ADGRB3", "ADGRL3", "ANKS1B", "ANO4", "ARNT2", "ASTN1", "ATP8A2", "CACNA1A", "CACNA1B", "CACNA1D", "CACNB4", "CADPS", "CELF4", "CNTN5", "CNTNAP2", "CNTNAP5", "CSMD1", "CSMD3", "DAB1", "DCLK1", "DLGAP1", "DLGAP2", "EPHA6", "ERC2", "FGF12", "FGF14", "FRMD4A", "FRMPD4", "GABBR2", "GABRB1", "GABRB2", "GABRG3", "GAD1", "GRIA1", "GRIA2", "GRID1", "GRIK2", "GRIN2A", "GRIN2B", "GRIP1", "GRM5", "HCN1", "IL1RAPL1", "KAZN", "KCNC2", "KCND2", "KCNJ3", "KCNQ3", "KCTD16", "KSR2", "LINC00632", "LRFN5", "LRRC7", "MAP2", "MDGA2", "MEG8", "MIAT", "MYT1L", "NMNAT2", "NRXN3", "NTRK2", "NXPH1", "OPCML", "PCLO", "PLCB1", "PTPRN2", "RAB3C", "RBFOX1", "RBFOX3", "RGS7", "RIMS1", "RIMS2", "RP11-384F7.2", "RP5-921G16.1", "RPS6KA2", "RTN1", "SCN1A", "SCN8A", "SLC35F1", "SLC44A5", "SLC4A10", "SLC6A1", "SNHG14", "SNRPN", "SNTG1", "SPOCK3", "SPTBN4", "STXBP5L", "SYN2", "SYT16", "TCF4", "TIAM1", "TMEM132B", "TMEM178B", "TRIM9", "UNC79", "XKR4", "ZNF385D", "ZNF536"))
gabaergic_interneuron
top20.idx<-to_plot_top20 %in% gabaergic_interneuron
gab20<-to_plot_top20[top20.idx]
gab20
DotPlot(all_sample, features = gab20, group.by = "tree.ident") + coord_flip()

lamp5_gabaergic_cortical_interneuron<-str_to_title(c("ADARB2", "SGCZ", "GRIN2A", "GRIK2", "FGF14", "PTPRT", "GRM5", "CSMD3", "ATP8A2", "GRIK1", "GRIP1", "SNHG14", "TMEM132B", "DAB1", "ZNF536", "LRRC7", "NXPH1", "TMEM132D", "GABRB2", "DLGAP1", "GRIN2B", "KAZN", "KCNC2", "MYT1L", "MTUS2", "GRIA4", "RBFOX1", "CNTNAP2", "GRIA1", "RGS7", "MYO16", "MGAT4C", "CACNA1B", "GABBR2", "SLC6A1", "CNTN5", "STXBP5L",
                                                     "SYN2", "IL1RAPL1", "FGF13", "KCNQ3", "KCTD16", "UNC5D", "XKR4", "CELF4", "SLC35F1", "GABRB1", "PLCB1", "ALK", "PTCHD4", "MDGA2", "FSTL5", "GRIA2", "DGKB", "CNTNAP4", "RBFOX3", "KCNIP1", "DNER", "NTRK2", "RP11-384F7.2", "RAB3C", "NRXN3", "RIMS1", "LSAMP", "ERBB4", "NRXN1", "LINGO2", "GRID1", "PTPRN2", "FRMPD4", "MCTP1", "SCN1A", "MACROD2", "LINC00632", "UNC79", "CACNA1A", "PRKCE", "ZMAT4", "SH3GL3", "DOCK3", "DOCK10", "GABRG3", "KIRREL3", "TMEM178B", "AC074363.1", "PDE8B", "ANKS1B", "RBMS3", "KCNJ3", "GRK3", "GPR158", "RIMBP2", "SYN3", "RAPGEF4", "CRACD", "CNTNAP5", "ADGRB3", "SOX2-OT", "GAD2", "GABRB3"))
lamp5_gabaergic_cortical_interneuron
top20.idx<-to_plot_top20 %in% lamp5_gabaergic_cortical_interneuron
lamp20<-to_plot_top20[top20.idx]
lamp20
DotPlot(all_sample, features = lamp20, group.by = "tree.ident") + coord_flip()

sncg_gabaergic_cortical_interneuron<-str_to_title(c("ADARB2", "GRIK2", "CNTN5", "MYT1L", "CSMD3", "DAB1", "RGS12", "SNTG1", "CNTNAP4", "ASIC2", "SNHG14", "CNTNAP2", "FSTL5", "GRIN2B", "CACNA1B", "NCAM2", "FRMPD4", "SHISA9", "NRXN1", "DSCAM", "STXBP5L", "FGF14", "GALNTL6", "CHRM3", "INPP4B",
                                                    "GABRB2", "LINGO2", "GRIK1", "GABBR2", "GABRB1", "NRXN3", "DLGAP1", "GRIA1", "ATP8A2", "KCNC2", "DGKB", "NPAS3", "PEX5L", "KCNQ5", "KCNQ3", "FRMD4A", "KCNJ3", "SPOCK3", "GRIA2", "SYNPR", "LRFN5", "LRRTM4", "SYN3", "EPHA6", "ERBB4", "ZNF536", "ZMAT4", "ADCY2", "PTPRN2", "PCSK2", "RBFOX3", "DOCK10", "KAZN", "HS6ST3", "MTUS2", "RIMS2", "OPCML", "KIAA1549L", "CELF4", "C8orf34", "NMNAT2", "GRM7", "SLC44A5", "PLCB1", "LINC00632", "MGAT4C", "PLD5", "KCND2", "CDH9", "CACNA1D", "RAB3C", "RP11-384F7.2", "DNER", "RPS6KA2", "PRELID2", "CHSY3", "SLC4A10", "ANKS1B", "ANO4", "ARPP21", "ARNT2", "DCC", "LRRC7", "RP11-...D23.1", "ZNF385D", "DOCK3", "TMEM178B", "ERC2", "CCNH", "SYN2", "L3MBTL4", "GABRB3", "GRIP1", "GRID1", "TNR"))
sncg_gabaergic_cortical_interneuron
top20.idx<-to_plot_top20 %in% sncg_gabaergic_cortical_interneuron
sncg20<-to_plot_top20[top20.idx]
sncg20
DotPlot(all_sample, features = sncg20, group.by = "tree.ident") + coord_flip()

vip_gabaergic_cortical_interneuron<-str_to_title(c("ADARB2", "CSMD1", "SNTG1", "SNHG14", "MYT1L", "SYNPR", "NRXN3", "GRM7", "GALNTL6", "DSCAM", "CHRM3", "ROBO2", "GRIK2", "LRRC7", "CNTNAP2", "NRXN1", "FGF14", "FRMD4A", "RGS12", "CSMD3", "ERBB4", "OPCML", "GRIN2B", "DLGAP1", "CACNA1B", "PTPRN2", "ATP8A2", "STXBP5L", "MDGA2", "KAZN", "PLCB1", "ANKS1B", "FRMPD4", "GABRB1", "ADGRB3", "DCLK1", "NRG3", "CNTN4", "RP11-384F7.2", "GABRB2", "ZNF536", "DAB1", "RGS7", "ROBO1", "PLD5", "PWRN1", "LINC00632", "GABBR2", "NLGN1", "GALNT13", "ANO4", "MEG8", "SLC44A5", "KCNJ3", "CACNA1D", "LRP1B", "GRIA2", "DLG2", "KCNH7", "L3MBTL4", "TMEM132B", "MTUS2", "SLC24A3", "RPS6KA2", "PLXNA4", "KCNMB2", "TRIM9", "LRRC4C", "SLC4A10", "GRIA1", "GRIP1", "MEG3", "NMNAT2", "IL1RAPL1", "FAM155A",
                                                   "TCF4", "ARNT2", "SYT16", "PAK3", "MAP2", "GABRG3", "RIMS2", "ASIC2", "KCNC2", "CDH10", "GRIN2A", "DOCK3", "CALN1","NBEA", "SEZ6L", "TAFA2", "DOCK10", "PPM1E", "DLX6-AS1", "UNC79", "FGF12", "NTRK3", "SNRPN", "CAMTA1", "KCNQ3"))
vip_gabaergic_cortical_interneuron
top20.idx<-to_plot_top20 %in% vip_gabaergic_cortical_interneuron
vip20<-to_plot_top20[top20.idx]
vip20
DotPlot(all_sample, features = vip20, group.by = "tree.ident") + coord_flip()

#glutamatergic<-str_to_title(c("CRACDL","RP11-213B3.1","IQCJ-SCHIP1","LINC02822","NWD2","SV2B","LINC00513","RP11-36B6.2","HSPA12A","MIR222HG","LY86-AS1","RP11-191L9.4","RP11-1263C18.2","IQSEC2","RP5-1015P16.1","AC114765.1","LINC01378","AC067956.1","NEFL","RP11-147G16.1","AC011288.2","TESPA1","RP11-170M17.1","COPG2IT1","CBLN2")) 
glutamatergic_neuron<-str_to_title(c("LRRC7","SATB2","MYT1L","ARPP21","NELL2","GRIN2B","NKAIN2","PTPRD","GRM5","KHDRBS2","CSMD1","MIR137HG","SLC4A10","RBFOX1","SNTG1","KCNH7","DLGAP1","GRM7","GRIA1","GABBR2","FRMPD4","FAT3","LRRTM4","ROBO2","SV2B","DSCAM","CSMD3","CACNA1B","TMEM132B","CAMK2A","ASIC2","CABP1","CACNA1E","CDH18","RALYL","GRIA3","DLGAP2","KCNH1","RBFOX3","PPP2R2C","SYT16","STXBP5L","CDH10","DCLK1","OLFM3","KIAA1549L","NETO1","GRIN1","PTPRN2","SPTBN4","SH3GL2","GABRB2","BASP1-AS1","HECW1","MMP16","GABRB3","VSNL1","LRFN5","CRACDL","NMNAT2","TNR","SCN8A","KCNC2","GABRG2","GRIA2","LINC01250","KCNQ3","CNTNAP5","DGCR5","KSR2","MDGA2","SYN2","SHISA9","SCN1A","SLIT1","DNM1","ANKS1B","AK5","LINC01122","FAM153CP","GABRG3","GRIN2A","KHDRBS3","RP1-232L24.3","ATP8A2","CHRM3","TAFA1","PEX5L","MIAT","KCNQ5","UNC5D","NDST3","CDH8","FRRS1L","SLC44A5","GPR158","LRFN2","OPCML","TAFA2","CACNG3"))
top20.idx<-to_plot_top20 %in% glutamatergic_neuron
glut20<-to_plot_top20[top20.idx]
glut20
DotPlot(all_sample, features = glut20, group.by = "tree.ident") + coord_flip()

near_projecting_glutamatergic_cortical_neuron<-str_to_title(c("ASIC2", "GRM5", "OLFM3", "SNTG1", "DLGAP2", "MDGA2", "CDH18", "RP11-419I17.1", "ZNF385D", "CSMD1", "GRIN2A", "TSHZ2", "ARPP21", "DCC", "PTPRN2", "HS3ST4", "CLSTN2", "GRIN2B", "RP11-...M17.1", "NKAIN2", "STXBP5L", "ITGA8", "FRMPD4", "SNHG14", "GRM3", "KHDRBS2", "CPNE4", "CSMD3", "GABRB2", "DCLK1", "GRM8", "SORCS2", "ERC2", "CHSY3", "CACNA1B", "TRHDE", "RP11-586K2.1", "MYT1L", "ADCY2", "EML6", "LRRC7", "SLC35F3", "SLC24A2", "NPSR1-AS1", "KHDRBS3", "RBFOX1", "CNTN4", "SLC4A10", "SPTBN4", "OPCML", "PAK5", "GABBR2", "TOX", "LRP1B", "CACNB4", "CADPS", "HS6ST3", "PWRN1", "SCN8A", "KSR2", "XKR6", "CDH8", "AK5", "GPR158", "PRR16", "DSCAM", "GABRG3", "HECW1", "LINC00632", "LDB2", "CACNA1A", "KCNQ3", "TLE4", "RGS7", "HTR2C", "GABRB1", "SH3RF3", "DAB1", "SV2B", "CHN1", "FAM135B", "SYN2", "PAK3", "CELF4", "LUZP2", "MEG8", "DPP10", "TRMT9B", "NAV3", "SH3GL2", "RP5-9...G16.1", "VWC2L", "FGF14", "MGAT4C", "NELL2", "MIAT", "PEX5L", "CNTNAP2", "KCNIP1", "KCNT2"))
top20.idx<-to_plot_top20 %in% near_projecting_glutamatergic_cortical_neuron
near20<-to_plot_top20[top20.idx]
near20
DotPlot(all_sample, features = near20, group.by = "tree.ident") + coord_flip()

Corticothalamic_projecting_glutamatergic_cortical_neuron<-str_to_title(c("ADAM23", "AK5", "ANKS1B", "ARPP21", "ASIC2", "ATP2B1", "CABP1", "CACNA1B", "CACNA1E", "CACNA2D3", "CACNB4", "CADPS", "CADPS2", "CDH18", "CELF4", "CHN1", "CLSTN2", "CNTNAP5", "CRACDL", "CSMD1", "DCLK1", "DLGAP2", "DNM1", "DPP10", "DSCAM", "ERC2", "FGF14", "FRMPD4", "GABBR2", "GABRB1", "GABRB2", "GABRG3", "GPR158", "GRIA2", "GRIA3", "GRIK4", "GRIN2A", "GRIN2B", "GRM3", "GRM5", "HECW1", "HS3ST4", "KALRN", "KCNH7", "KCNQ3", "KCNQ5", "KHDRBS2", "KHDRBS3", "KIAA1549L", "KSR2", "LDB2", "LINC00632", "LINC01122", "LRFN5", "LRP1B", "LRRC7", "LRRTM4", "MCTP1", "MDGA2", "MMP16", "MYT1L", "NELL2", "NKAIN2", "NRG1", "OLFM3", "OPCML", "PAK5", "PCLO", "PDZRN4", "PEX5L", "PRKCB", "PTPRD", "PTPRN2", "RALYL", "RAPGEF4", "RBFOX1", "RBFOX3", "RYR3", "SCN8A", "SEMA3E", "SH3GL2", "SH3GL3", "SHISA9", "SLC35F1", "SLC4A10", "SLIT1", "SNHG14", "SNTG1", "SORCS1", "SORCS3", "SPTBN4", "STXBP5L", "SV2B", "SYN2", "TLE4", "TMEFF2", "TMEM132B", "UNC79", "UNC80", "VSNL1"))
top20.idx<-to_plot_top20 %in% Corticothalamic_projecting_glutamatergic_cortical_neuron
cort20<-to_plot_top20[top20.idx]
cort20
DotPlot(all_sample, features = cort20, group.by = "tree.ident") + coord_flip()

l6b_glutamatergic_cortical_neuron<-str_to_title(c("AK5", "ARPP21", "ASIC2", "BASP1-AS1", "BRINP2", "CABP1", "CACNA1B", "CACNA1E", "CACNA2D3", "CACNB4", "CADPS", "CDH10", "CDH18", "CELF2", "CELF4", "CHN1", "CHRM3", "CNTNAP5", "CSMD1", "DCLK1", "DLGAP1", "DLGAP2", "DNM1", "DPP10", "DSCAM", "ERC2", "FGF14", "FOCAD", "FRMPD4", "FUT9", "GABBR2", "GABRB2", "GABRB3", "GABRG3", "GARNL3", "GAS7", "GPR158", "GRIA1", "GRIA2", "GRIA3", "GRIK2", "GRIN2B", "GRM5", "HECW1", "HS3ST4", "KALRN", "KCNJ6", "KHDRBS2", "KHDRBS3", "KIAA1217", "KIAA1549L", "KIFAP3", "KSR2", "LINC00632", "LINC01122", "LMO3", "LRP1B", "LRRC7", "MCTP1", "MDGA2", "MIAT", "MIR137HG", "MMP16", "MYT1L", "NELL2", "NKAIN2", "NMNAT2", "NOS1AP", "OLFM3", "PAK5", "PCDH11X", "PCSK5", "PDZRN4", "PPP2R2C", "PTPRD", "PTPRN2", "PTPRR", "RALYL", "RBFOX1", "RBFOX3", "RGS7", "ROBO2", "SATB2", "SCN8A", "SEZ6L", "SH3GL2", "SH3GL3", "SLC24A2", "SLC35F1", "SLC4A10", "SNHG14", "SORCS3", "SPTBN4", "STXBP5-AS1", "STXBP5L", "SV2B", "SYT16", "TLE4", "TMEM132B", "XKR4"))
top20.idx<-to_plot_top20 %in% l6b_glutamatergic_cortical_neuron
l6b20<-to_plot_top20[top20.idx]
l6b20
DotPlot(all_sample, features = l6b20, group.by = "tree.ident") + coord_flip()

astrocyte<-str_to_title(c("APOE", "AQP1", "CRYAB", "EDNRB", "GFAP", "GPC5", "HPSE2", "ITPRID1", "LGR6", "LINC00499", "LINC00836", "OBI1-AS1", "RP11-...K15.2", "RP11-297I23.1", "RP11-517I3.1", "SLC14A1", "TNC"))
top20.idx<-to_plot_top20 %in% astrocyte
as20<-to_plot_top20[top20.idx]
as20
DotPlot(all_sample, features = astrocyte, group.by = "tree.ident") + coord_flip()

microglial_cell<-str_to_title(c("AIF1", "APBB1IP", "CCDC26", "CX3CR1", "DOCK8", "F13A1", "HIF1A-AS3", "ITGAX", "LINC00278", "LINC02642", "LINC02712", "LNCAROD", "MRC1", "OLR1", "P2RY12", "SCIN", "SRGN", "TTR"))
top20.idx<-to_plot_top20 %in% microglial_cell
mic20<-to_plot_top20[top20.idx]
mic20
DotPlot(all_sample, features = mic20, group.by = "tree.ident") + coord_flip()

cerebral_cortex_endothelial_cell<-str_to_title(c("ABCB1", "ABCG2", "ADGRF5", "ADGRL4", "ADIPOR2", "AGFG1", "ANKS1A", "ANO2", "APOLD1", "ARHGAP29", "ARHGAP31", "ARL15", "ATP10A", "BMPR2", "BTNL9", "CADPS2", "CCNY", "CDYL2", "CEMIP2", "CHSY1", "CLDN5", "CMTM8", "CNOT6L", "CPNE8", "CRIM1", "DOCK1", "DOCK9", "EGFL7", "ELOVL7", "EPAS1", "EPB41L4A", "ERG", "ESYT2", "FLI1", "FLT1", "GALNT15", "GALNT18", "GRB10", "HERC2", "HIF1A-AS3", "HIPK3", "IGF1R", "IL4R", "IRAK3", "ITGA1", "LDLRAD3", "LEF1", "LHFPL6", "LINC00472", "LMBR1", "LRCH1", "MCF2L", "MECOM", "MEF2A", "MRTFB", "MTUS1", "MYRIP", "NEDD9", "NOSTRIN", "NXN", "PECAM1", "PICALM", "PIK3R3", "PLEKHG1", "PLXNA2", "PODXL", "PON2", "PPFIBP1", "PREX2", "PRKCH", "PTPRB", "PTPRG", "RAPGEF1", "RAPGEF2", "RBMS2", "RNF144B", "RPGR", "RUNDC3B", "SCARB1", "SEC14L1", "SGPP2", "SLC16A1", "SLC1A1", "SLC2A1", "SLC39A10", "SLC7A1", "SLC7A5", "SLCO2B1", "SLCO4A1", "SORT1", "ST6GAL1", "ST6GALNAC3", "ST8SIA6", "TACC1", "TBC1D4", "TGM2", "THSD4", "USP6NL", "VWF", "WWTR1"))
top20.idx<-to_plot_top20 %in% cerebral_cortex_endothelial_cell
cer20<-to_plot_top20[top20.idx]
cer20
DotPlot(all_sample, features = cer20, group.by = "tree.ident") + coord_flip()

central_nervous_system_macrophage<-str_to_title(c("ABCC4", "AC008697.1", "ADAM28", "ANKRD44", "AOAH", "APBB1IP", "ARHGAP15", "ARHGAP22", "ARHGAP24", "ATP8B4", "BMP2K", "CHST11", "CSF1R", "CSGALNACT1", "CX3CR1", "CYFIP1", "DISC1", "DLEU1", "DOCK10", "DOCK2", "DOCK4", "DOCK8", "ELMO1", "EPB41L2", "FRMD4A", "HS3ST4", "INPP5D", "ITPR2", "KCNQ3", "LDLRAD4", "LHFPL2", "LINC01374", "LINC02232", "LINC02798", "LNCAROD", "LPAR6", "LPCAT2", "LRMDA", "LRRK1", "MAML3", "MEF2A", "MEF2C", "MGAT4A", "P2RY12", "PALD1", "PLXDC2", "PREX1", "RASGEF1C", "RP11-358F13.1", "RUNX1", "SFMBT2", "SLC1A3", "SLC9A9", "SLCO2B1", "SORL1", "SRGAP2", "SRGAP2B", "ST6GAL1", "ST6GALNAC3", "SYNDIG1", "TBXAS1"))
top20.idx<-to_plot_top20 %in% central_nervous_system_macrophage
cen20<-to_plot_top20[top20.idx]
cen20
DotPlot(all_sample, features = cen20, group.by = "tree.ident") + coord_flip()

SaveObject(all_sample, "seurat_obj_after_PCA")

all_sample <- ReadObject("seurat_obj_after_PCA")


#DotPlot(all_sample, features = to_plot, group.by = "tree.ident") + coord_flip()

#markers<-c("Angpt4", "Arhgap15", "Mobp","Dock6", "Pdgfra", "Pdgfrb", "Eng", "Syt6","Slc9c1", "Mbp")
# syt6 is glutamatergic from cellxgene canonical for glut nuerons.
#DotPlot(all_sample, features=markers, group.by="tree.ident")+coord_flip()

#DimPlot(all_sample, reduction = "umap", label = TRUE) + NoAxes() + NoLegend()

new_ids <- c("Cer.Cort.Endot.Cell", "Slc7a10", "Adgrl4", "Sst", "Cort.Proj.Glut.Cort. + L6b.Glut.Cort.Neuron", "Near.Proj.Glut.Cort.Neuron","Cftr","Opalin","Microglial.Cell", "LOC102555705",
             "Ntf3",  "LOC103692065", "Lamp5.Gab.Cort +Sncg.Gab.Cort+Vip.Gab.Cort.Interneuron")

new_id_list <- list("Cer.Cort.Endot.Cell" = 1, "Slc7a10" = 2, "Adgrl4" = 3, "Sst" = 4,"Cort.Proj.Glut.Cort. + L6b.Glut.Cort.Neuron" = 5, "Near.Proj.Glut.Cort.Neuron" = 6,"Cftr"=7, "Opalin" = 8,"Microglial.Cell" = 9, "LOC102555705" = 10,"Ntf3"=11,
                    "LOC103692065" = 12, "Lamp5.Gab.Cort +Sncg.Gab.Cort+Vip.Gab.Cort.Interneuron"=c(13,14))

#new_ids <- c("Astrocyte", "Loc103692025", "Reln", "Near.Proj.Glut.Cort.Neuron","Plp1","Il1rapl2","Loc108353456",
#             "Glut.Neuron", "Cort.Proj.Glut.Cort. + L6b.Glut.Cort.Neuron", "Cer.Cort.Endot.Cell", "Microglial.Cell+Cent.Nerv.Sys.Macrophage", 
#             "Gab+Lamp5.Gab.Cort +Sncg.Gab.Cort+Vip.Gab.Cort.Interneuron")

#new_id_list <- list(Astrocyte = 1, Loc103692025 = 2, Reln = 4,Near.Proj.Glut.Cort.Neuron = 6,Plp1 = 8,Il1rapl2 = 9,Loc108353456=11,
#                    Glut.Neuron = c(7,12), Cort.Proj.Glut.Cort_L6b.Glut.Cort.Neuron = 5, Cer.Cort.Endot.Cell = 3, Microglial.Cell_Cent.Nerv.Sys.Macrophage = 10, 
#                    Gab_Lamp5.Gab.Cort_Sncg.Gab.Cort_Vip.Gab.Cort.Interneuron=c(13,14))
#new_id_list[[8]]
#which(is.na(all_sample@meta.data$tree.ident))
all_sample@meta.data$collapsed<-NA
for (i in 1:length(new_id_list)) {
    ind <- which(all_sample@meta.data$tree.ident %in% new_id_list[[i]])
    all_sample@meta.data$collapsed[ind] <- names(new_id_list)[i]

  }

#which(is.na(all_sample@meta.data$collapsed))

all_sample@meta.data$collapsed <- factor(
  all_sample@meta.data$collapsed, levels = names(new_id_list), ordered = TRUE)
Idents(all_sample) <- all_sample@meta.data$collapsed

names(new_ids) <- levels(all_sample)
all_sample <- RenameIdents(all_sample, new_ids)
pdf (file="umap_mclover3.pdf", width=6, height=6)
DimPlot(all_sample, reduction = "umap", label = TRUE) + NoLegend() #NoAxes() + NoLegend()
dev.off()
