
args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0) {
	stop("No arguments supplied\n", call.=FALSE)
}

library(Seurat)
library(domino2)
library(dplyr)

working_dir <- "results"
grouping <- as.character(args[1])
group <- as.character(args[2])
cluster.col <- as.character(args[3])

seu <- readRDS("data/3_ser.RDS")
seu$Patient_Tissue <- paste0(seu$Patient, "_", seu$Tissue)
seu$SnC <- ifelse(is.na(seu$projectR0.01), "not_SnC", "SnC")
seu$cell_type_SnC <- seu$cell_type
seu$cell_type_SnC[seu$cell_type_SnC %in% c("Fibroblast", "Pericyte_SMC", "SMC", "Lymphatic_endo", "Endo_vein2")] <- paste0(seu$cell_type_SnC[seu$cell_type_SnC %in% c("Fibroblast", "Pericyte_SMC", "SMC", "Lymphatic_endo", "Endo_vein2")], "_", seu$SnC[seu$cell_type_SnC %in% c("Fibroblast", "Pericyte_SMC", "SMC", "Lymphatic_endo", "Endo_vein2")])
seu <- seu[,! seu$cell_type %in% c("Damaged_endo", "Damaged_Epi", "Damaged_mye", "Doublet")]



if(!dir.exists(file.path(working_dir, "domino", cluster.col))) {
        dir.create(file.path(working_dir, "domino", cluster.col))
}

if(!dir.exists(file.path(working_dir, "domino", cluster.col, group))) {
	dir.create(file.path(working_dir, "domino", cluster.col, group))
}

cat("Running domino2 on ", grouping, " ", group, "with ",  cluster.col, " clustering\n")

#pyscenic results
auc <- read.table(file.path(working_dir, "scenic", "fibroid_auc.csv"), header = T, row.names = 1, stringsAsFactors = F, sep = ",")
colnames(auc) <- sapply(colnames(auc), gsub, pattern = "\\.\\.\\.", replacement = "")
regulon_df <- read.csv(file.path(working_dir, "scenic", "fibroid_regulons.csv"), skip = 2)
regulon_colnames <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity", "Annotation", "Context", "TargetGenes", "RankAtMax")
colnames(regulon_df) <- regulon_colnames
regulon_ls <- create_regulon_list_scenic(regulon_df)

#seurat data
seu.tx <- seu[,seu[[grouping]] == group]
seu.tx
counts <- GetAssayData(seu.tx, assay = "RNA", slot = "counts")
dim(counts)
z_scores <- GetAssayData(seu.tx, assay = "RNA", slot = "scale.data")
dim(z_scores)
clusters <- seu.tx[[cluster.col]]
str(clusters)
clusters <- as.factor(clusters[,1])
names(clusters) <- colnames(seu.tx)
table(clusters)
cells <- colnames(seu.tx)
length(cells)

auc <- auc[cells, ]

cellphonedb <- readRDS(file.path(working_dir, "domino", "cpdb_rl_map.rds"))
head(cellphonedb)

dom <- create_domino(
  rl_map = cellphonedb,
  features = t(auc),
  counts = counts,
  z_scores = z_scores,
  clusters = clusters,
  tf_targets = regulon_ls,
  use_clusters = T,
  use_complexes = T,
  remove_rec_dropout = F
) %>% build_domino(
  dom = .,
  min_tf_pval = .001,
  max_tf_per_clust = Inf,
  max_rec_per_tf = Inf,
  rec_tf_cor_threshold = .25
)

saveRDS(dom, file.path(working_dir, "domino", cluster.col, group, paste0(group, ".domino.rds")))

rm(list = ls())
warnings()
sessionInfo()


