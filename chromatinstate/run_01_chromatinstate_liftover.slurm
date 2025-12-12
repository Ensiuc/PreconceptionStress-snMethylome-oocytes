#!/bin/bash
#SBATCH --job-name=chromatinstate
#SBATCH --output=jupyter_chrmatin_output.log
#SBATCH --error=jupyter_chromatin_error.log
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G

#module load python  # Load Python if needed
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ~/miniconda3/envs/allcools  # Activate your environment if needed

#python 01-chromatinestate-liftingregions.py  # Run the converted Python script

python 01-lifthyperhypoonly.py 
