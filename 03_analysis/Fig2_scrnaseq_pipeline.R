library(Seurat)
library(DoubletFinder)
library(harmony)
library(RColorBrewer)
library(dplyr)

result_path <- "output/"

sample <- c("Tumor2-23", "Tumor2-25-1", "Tumor2-25-2", "Tumor2-27-1", "Tumor2-27-2", "Tumor3-1-1",
            "Tumor3-1-2", "Tumor3-3-1", "Tumor3-3-2", "Tumor3-12-1", "Tumor3-12-2", "Label_NC_1", "Label_NC_2")

cell_markers <- list(
  Macrophage  = c("Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "C1qa", "C1qb", "Apoe", "Cd68"),
  cancer_cell = c("Epcam", "Krt8", "Krt18", "Vim"),
  PCT         = c("Lrp2", "Slc5a12", "Slc13a3", "Slc16a9", "Inmt"),
  ENDO        = c("Plat", "Emcn", "Plpp1", "Ehd3", "Nrp1", "Kdr"),
  NK_T        = c("Ccl5", "Gzma", "Nkg7", "Gzmb", "Ltb", "Cxcr6", "Il7r"),
  PST         = c("Aadat", "Kap", "Napsa", "Slc22a13"),
  Neutrophil  = c("S100a8", "S100a9", "Retnlg", "Il1b", "Ngp"),
  TAL         = c("Slc12a1", "Umod", "Egf", "Wfdc15b", "Sostdc1"),
  Fibroblast  = c("Col1a1", "Col1a2", "Dcn", "Fap", "Pdpn"),
  CDIC_CDPC   = c("Aqp2", "Hsd11b2", "Fxyd4", "Atp6v1g3", "Atp6v0d2", "Slc26a4")
)

rm_doublet <- function(name = NULL, input = NULL, dim.usage = 30) {
  colnames(input) <- paste0(name, "_", colnames(input))
  seurat_obj <- CreateSeuratObject(input, project = name, min.cells = 3, min.features = 50)
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^mt-", col.name = "percent.mt")
  seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:dim.usage)

  Find_doublet <- function(data) {
    sweep.res.list <- paramSweep(data, PCs = 1:dim.usage, sct = TRUE)
    sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn          <- find.pK(sweep.stats)
    nExp_poi       <- round(0.05 * ncol(data))
    p              <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK))
    data           <- doubletFinder(data, PCs = 1:dim.usage, pN = 0.25, pK = p,
                                    nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    colnames(data@meta.data)[ncol(data@meta.data)] <- "doublet_info"
    return(data)
  }

  seurat_obj <- Find_doublet(seurat_obj)
  seurat_obj <- subset(seurat_obj, subset = doublet_info == "Singlet")
  seurat_obj@meta.data$library <- name

  c_idx <- grep("pANN_", colnames(seurat_obj@meta.data))
  seurat_obj@meta.data <- seurat_obj@meta.data[, -c_idx]

  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

  pdf(paste0(result_path,  name, "_before_quality.pdf"))
  print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, group.by = "orig.ident"))
  dev.off()

  max_mt_value      <- quantile(seurat_obj$percent.mt,    probs = 0.975)
  min_Feature_value <- quantile(seurat_obj$nFeature_RNA,  probs = 0.025)
  max_Feature_value <- quantile(seurat_obj$nFeature_RNA,  probs = 0.975)

  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA > min_Feature_value &
                                nFeature_RNA < max_Feature_value &
                                percent.mt   < max_mt_value)

  pdf(paste0(result_path,  name, "_after_quality.pdf"))
  print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, group.by = "orig.ident"))
  dev.off()

  saveRDS(seurat_obj, paste0(result_path,  name, ".rds"))
  return(seurat_obj)
}

rds_list <- list()
for (sam in sample) {
  rds_list[[sam]] <- readRDS(paste0(result_path,  sam, ".rds"))
}

features    <- SelectIntegrationFeatures(object.list = rds_list, nfeatures = 3000)
rds_list    <- PrepSCTIntegration(object.list = rds_list, anchor.features = features)
all.anchors <- FindIntegrationAnchors(object.list = rds_list, normalization.method = "SCT",
                                      anchor.features = features)
seurat_obj  <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT")

seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "library")
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:20,
                      reduction.name = "umap", min.dist = 0.03, n.neighbors = 5L)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)

options(repr.plot.height = 5, repr.plot.width = 12)
VlnPlot(object = seurat_obj, features = c("nCount_RNA", "nFeature_RNA"),
        group.by = "seurat_clusters", pt.size = 0.1)

options(repr.plot.height = 5, repr.plot.width = 5)
custom_color <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_obj@meta.data$seurat_clusters)))
print(DimPlot(seurat_obj, cols = custom_color, label = TRUE))

new.cluster.ids <- c("Macrophage", "Macrophage", "Macrophage", "cancer cell", "PCT", "cancer cell",
                     "cancer cell", "PCT", "PCT", "Macrophage", "ENDO",
                     "NK_T", "PST", "cancer cell", "Macrophage", "cancer cell",
                     "Neutrophil", "PCT", "TAL", "Macrophage", "Fibroblast",
                     "cancer cell", "CDIC_CDPC", "Macrophage", "cancer cell")

Idents(seurat_obj) <- seurat_obj$seurat_clusters
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$celltype1 <- Idents(seurat_obj)

custom_color <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_obj@meta.data$celltype1)))
DimPlot(seurat_obj, reduction = "umap", label = TRUE, cols = custom_color,
        pt.size = 0.3, group.by = "celltype1")


###
library(Seurat)
library(future)
library(copykat)

plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^5)
set.seed(123)

result_path <- "copykat_result/"
ALL <- readRDS("1_main_seu.rds")

Idents(ALL) <- ALL$sample_name
sample_index <- c("Tumor2-23", "Tumor2-25-1", "Tumor2-25-2", "Tumor2-27-1", "Tumor2-27-2", "Tumor3-1-1",
            "Tumor3-1-2", "Tumor3-3-1", "Tumor3-3-2", "Tumor3-12-1", "Tumor3-12-2", "Label_NC_1", "Label_NC_2")

for (sample in sample_index){
  print(sample)
  subcells <- WhichCells(ALL, idents = sample)
  copykat_res <- subset(ALL, cells = subcells)
  print(table(copykat_res@meta.data$sample_name))
  copykat_res_rawdata <- as.data.frame(copykat_res[['RNA']]$counts)
  copykat.test <- copykat(rawmat=copykat_res_rawdata, 
                          id.type="S", 
                          cell.line="no", 
                          ngene.chr=5, 
                          win.size=25, 
                          KS.cut=0.15, 
                          sam.name=sample, 
                          distance="euclidean", 
                          n.cores=50,genome="mm10")
}
