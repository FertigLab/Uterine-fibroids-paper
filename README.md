# Senescence and dominoSignal analysis for the paper "Senescent cell networks link matrix remodeling and vascular dysfunction in human uterine fibroids"

### Data and preprocessing
The single-cell RNA sequencing raw data used in this analysis will be available for download on NCBI GEO [upload pending]. The preprocessing and annotation for the scRNA-seq was performed as described in [fibroid_reproducibility](Elisseeff-Lab/fibroid_reproducibility: Code to reproduce various fibroid single cell analyses). 

### To make a local copy of this repository

```sh
git clone https://github.com/FertigLab/Uterine-fibroids-paper.git
```

### Install dominoSignal

Download v0.2.2 from [dominoSignal releases](https://github.com/FertigLab/dominoSignal/releases)
```R
install.packages("domino2-0.2.2-alpha.tar.gz", repos = NULL, type = "source")

#OR

remotes::install_github("FertigLab/dominoSignal@v0.2.2-alpha")
```

### Analysis
1. Transfer learning with projectR to score the dataset for the SenSig.
```sh
├── 1_projectR
│   ├── UMAP_celltype_tissue.Rmd - Fig 1F
│   ├── transfer_learning.rmd
│   ├── SnC_plots.Rmd - Fig 1H 
│   ├── plot_SnC_by_patient_tissue.rmd
│   ├── SnC_related_genes.umap.Rmd
│   └── projection_drivers.Rmd
```
Generates SenSig projection weights and projection p-values for all the cells in the dataset and makes plots. Also generates a list of genes in the signature that distinguish the senescent and non-senescent populations as projection drivers.

2. Analysis of the senescent populations.
```sh
├── 2_SnC_cell_analysis
│   ├── SnC_cell_analysis.Rmd
│   ├── SnC_cell_analysis_volcano_barplot.Rmd 
│   ├── heatmap_SnC_vs_nonSnC_FC.Rmd - Fig 2A
│   ├── heatmap_SnC_vs_nonSnC_FC.secreted.Rmd - Fig 2B
│   └── Summary_GSEA_plot.Rmd - Fig 2C
```
Generates a list of genes that are differentially expressed in the senescent cells as compared to the non-senescent cells by cluster. Also adds CellphoneDB(v4, included in ref/cpdb/v4.0.0) scereted genes, MatrixDB, surfaceome annotations. Includes scripts to generate plots in Figure 2.

3. Cell-cell communication analysis with dominoSignal
```sh
├── 3_domino
│   ├── fibroid_scenic.sh
│   ├── create_RL_map.Rmd
│   ├── domino.R
│   ├── submit_domino.sh
│   ├── domino_list.R
│   ├── differential_linkages.rec_lig_activation.Rmd
│   ├── SnC_vascular_to_stromal_vascular_network.Rmd
│   ├── SnC_stromal_to_immune.Rmd
│   ├── differential_TF_activation.SnC_vs_nonSnC.Rmd
│   ├── Mechano_TFs_plot.Rmd
│   ├── plot_TF_targets.Rmd
│   ├── all_signaling_to_SnC.receptors_connected_to_SRF_TEAD4.Rmd
│   ├── integrins_only_without_complexes
│   │   ├── domino.R
│   │   ├── submit_domino.sh
│   │   ├── domino_list.R
│   │   ├── differential_linkages.rec_lig_activation.Rmd
│   │   └── collagen_network_from_Fibroblasts_Pericytes.Rmd
│   └── TF_activity_on_umap.Rmd
```

Run [SCENIC](https://doi.org/10.1038/nmeth.4463) analysis on the dataset
```sh
sbatch fibroid_scenic.sh
```

Create a receptor-ligand interaction dataframe and launch dominoSignal by tissue/patient_tissue using domino.R as below:
```sh
create_RL_map.Rmd 
sbatch submit_domino.sh "SAMPLEID_METADATA_COLUMN_NAME" "SAMPLEID" cell_type_SnC
example: sbatch submit_domino.sh Patient_Tissue sampleID_Fibroid cell_type_SnC
```
Generates dominoSignal analysis with the inter- and intra- cellular interactions at a cluster level encoded in an object.

For the dominoSignal analyses performed by patient_tissue, create a [linkage summary object](https://academic.oup.com/bioinformatics/article/42/3/btag089/8499660) and run differential receptor-ligand signaling comparing fibroid and myometrium samples.
```sh
domino_list.R
differential_linkages.rec_lig_activation.Rmd
```
This returns a differential signaling table with a Fisher's exact test result and tabulation of the number of samples a linkage is observed in per group. 

Explore signaling networks highlighted in the paper.
```sh
SnC_vascular_to_stromal_vascular_network.Rmd - Fig 3 I,J
SnC_stromal_to_immune.Rmd - Fig 4 E,G
```

Compare transcription factor activation levels (using the TF AUC values from SCENIC analysis) between senescent and non-senescent cells.
```sh
differential_TF_activation.SnC_vs_nonSnC.Rmd
```

Explore Mechanosensing TFs.
```sh
Mechano_TFs_plot.Rmd - Fig 5D
plot_TF_targets.Rmd - Fig 5 E,F
all_signaling_to_SnC.receptors_connected_to_SRF_TEAD4.Rmd - Fig 5G
```

Run dominoSignal without complexes and explore signaling using scripts under integrins_only_without_complexes. Usage and description similar to above.  
```sh
collagen_network_from_Fibroblasts_Pericytes.Rmd - Fig 5H
```
