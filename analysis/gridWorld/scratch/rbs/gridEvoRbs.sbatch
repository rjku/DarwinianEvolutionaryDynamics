#!/bin/sh

#SBATCH --output=rbsGridEvo-%j-slurmLog.out

#SBATCH --time=12:00:00
#SBATCH --ntasks=200

module load anaconda3
module load julia

folderName="/scratch/riccardorao/gridWorld/rbs/"
jobTag="rbsGridEvo"
julia generateRbsTraj.jl $SLURM_NTASKS $folderName $SLURM_JOB_ID $jobTag

# ############################

# #SBATCH --time=10:00:00
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node=20

# #SBATCH --cpus-per-task=1
# #SBATCH --exclusive
# #SBATCH --exclusive=user

# #SBATCH --output=%j-SlurmLog.out
