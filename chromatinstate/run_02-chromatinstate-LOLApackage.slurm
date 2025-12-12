#!/bin/bash
#SBATCH --job-name=chromatinstate
#SBATCH --output=rbase_chromatin_output.log
#SBATCH --error=rbase_chromatin_error.log
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G


###################
# Run Information #
###################

start=`date +%s`

hostname

THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_CPU} / 1024)

echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

####################
# Load Environment #
####################

# Activate the conda environment
conda activate /quobyte/lasallegrp/programs/.conda/envs/DMRichR_R4.2


# Run with defaults

call="Rscript \
--vanilla \
/quobyte/lasallegrp/Ensi/project/oocytes-macaque/allcools/allcools-subsampled/alloocytes/results-pooled/chromatinstate/02-chromatinstate-LOLApackage.R"

echo $call
eval $call

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
