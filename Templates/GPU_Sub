#!/bin/bash
#SBATCH --job-name="{Sim_Name}"
#SBATCH --output="lammpsgpu.%j.%N.out"
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --no-requeue
#SBATCH --gres=gpu:4
#SBATCH -t 048:00:00
#SBATCH -A csd467
cd {path}

module unload mvapich2_ib
module unload intel
module load intel/2015.2.164
module load mvapich2_ib
module load cuda

ibrun -np 8 /share/apps/gpu/lammps/lmp_cuda_mpi -sf gpu -pk gpu 4 -in in.{Sim_Name} -log log.{Sim_Name}
