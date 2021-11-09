#!/bin/bash
#SBATCH -N 12
#SBATCH --qos=regular
#SBATCH --time=02:00:00
#SBATCH --constraint=haswell
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=kbrauer@lbl.gov
#SBATCH --mail-type=ALL

export MPICH_MAX_THREAD_SAFETY=multiple # NEED THIS OR WILL FAIL

# cd $SLURM_SUBMIT_DIR
# export OMP_PLACES=cores
# export OMP_PROC_BIND=true
# export OMP_NUM_THREADS=4

# srun --cpu_bind=cores ./Castro2d.gnu.haswell.DEBUG.MPI.ex inputs.2d.sph2
srun --cpu_bind=cores ./Castro3d.gnu.haswell.DEBUG.MPI.ex inputs.3d.sph
