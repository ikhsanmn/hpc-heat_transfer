#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
cd $PWD
mpirun -np 16 --oversubscribe python3 py_LBM3D_coupled_without_cycle.py

