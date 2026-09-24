#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G

module load conda_R/4.3
module list

R CMD BATCH "--args $1 $2 $3" domino.R ${3}/${2}/${2}.domino.Rout
