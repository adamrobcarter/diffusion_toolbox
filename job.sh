#!/bin/bash
#SBATCH --output=slurm%j.out
#SBATCH --job-name=toolbox
#SBATCH --mail-user=adam.carter@sorbonne-universite.fr
#SBATCH --mail-type=ALL
#SBATCH --partition=std
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=64000M
set -e # stop on error
#export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK}}
## Load software modules
# module load cuda/12.8
## Move data (if necessary) and launch application
## (use 'srun' for launching mpi applications)
cd ~/toolbox

python -m isf.calc_f ld_hydro_dpstokes_0.114_L5120_t1h_64

## Finish gracefully
exit 0