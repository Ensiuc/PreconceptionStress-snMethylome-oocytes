#!/bin/bash
#
#SBATCH --job-name=DM.R
#SBATCH --ntasks=4 # Number of cores/threads
#SBATCH --mem=300G # Ram in Mb
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --time=5-00:00:00

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
##########################################################################################

###################
# Run Information #
###################

start=`date +%s`

hostname

THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_CPU} / 1024)

echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

################
# Load Modules #
################

module load R/4.1.0
module load homer

#source ~/.bashrc
#conda activate /quobyte/lasallegrp/programs/.conda/envs/DMRichR4.2_Updated/
#conda activate /quobyte/lasallegrp/programs/.conda/envs/DMRichR_R4.2/



########
# DM.R #
########

call="Rscript \
--vanilla \
--{pathto}/DM_R4.2.R \
--genome rheMac10 \
--coverage 1 \
--perGroup '0.75' \
--minCpGs 5 \
--maxPerms 10 \
--maxBlockPerms 10 \
--cutoff '0.05' \
--testCovariate Group \
--adjustCovariate 'Year' \
--sexCheck FALSE \
--GOfuncR TRUE \
--EnsDb FALSE \
--cores 2"

echo $call
eval $call

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
