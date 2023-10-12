#!/bin/bash

#SBATCH -A teska
#SBATCH -J NSForest
#SBATCH --mem 100000M
#SBATCH -D /nfs/home/students/teska/omnideconv_input/work/nsforest

srun ~/masterthesis/scripts/nextflow/nsf_submission.sh
