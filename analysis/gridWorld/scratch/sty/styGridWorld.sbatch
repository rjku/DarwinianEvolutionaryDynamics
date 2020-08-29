#!/bin/sh

#SBATCH --output=styGridEvo-%j-slurmLog.out

#SBATCH --time=24:00:00
#SBATCH --ntasks=1000

module load anaconda3
module load julia

folderName="/scratch/riccardorao/gridWorld/sty/"
jobTag="styGridEvo"
julia generateStyTraj.jl $SLURM_NTASKS $folderName $SLURM_JOB_ID $jobTag

# ############################

# #SBATCH --time=10:00:00
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node=20

# #SBATCH --cpus-per-task=1
# #SBATCH --exclusive
# #SBATCH --exclusive=user

# #SBATCH --output=%j-SlurmLog.out
