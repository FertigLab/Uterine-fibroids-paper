#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G

module load r/4.4.2
module list

RSCRIPT="domino.R"
ROUT="results/domino/${3}/${2}/${2}.domino.Rout"

Rscript --verbose --no-save $RSCRIPT $1 $2 $3 > $ROUT

exit 0
