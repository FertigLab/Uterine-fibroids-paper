#!/bin/bash

#SBATCH --cpus-per-task 32
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=384G

module load singularity

# Change to current directory:
#cd run_scenic/fibroid/

# Running SCENIC on combined fibroid data
SIF_DIR="./" # Singularity image directory
SCENIC_DIR="./" # For human data
PREFIX="fibroid" # Project name or identifier
INPUT_LOOM="fibroid_v3_counts.loom" # Loom file with counts
OUTPUT_DIR="./" # Directory for outputs
NUM_CORE=31 # Number of cores to use

# Co-expression modules are used to quantify gene-transciption factor adjacencies.
singularity exec "aertslab-pyscenic-0.12.1.sif" \
    pyscenic grn \
    "${INPUT_LOOM}" "allTFs_hg38.txt" --output "${PREFIX}_adj.tsv" \
    --num_workers ${NUM_CORE} --seed 42

# ctx: construct TF regulons with pruning based on TF motiff enrichment
singularity exec "aertslab-pyscenic-0.12.1.sif" pyscenic ctx \
    "${PREFIX}_adj.tsv" "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    --annotations_fname "motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
    --expression_mtx_fname "${INPUT_LOOM}" --mode "dask_multiprocessing" --output "${PREFIX}_regulons.csv" \
    --num_workers ${NUM_CORE}

#### aucell: calculate TF activity scores
# AUCell calculations
singularity exec "aertslab-pyscenic-0.12.1.sif" pyscenic aucell \
    "${INPUT_LOOM}" "${PREFIX}_regulons.csv" -o "${PREFIX}_auc.csv"