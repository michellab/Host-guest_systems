#!/bin/bash
#SBATCH --job-name=pmd_cG4_van
#SBATCH -o somd-array-gpu-%A.%a.out
#SBATCH -p gpu
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --time 24:00:00
#SBATCH --array=0-10
module load cuda/9.2
echo "CUDA DEVICES:" $CUDA_VISIBLE_DEVICES
lamvals=(0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 )
lam=${lamvals[SLURM_ARRAY_TASK_ID]}
echo "lambda is: " $lam
mkdir lambda-$lam
cd lambda-$lam
export OPENMM_PLUGIN_DIR=/export/users/sofia/sire-jm.app/lib/plugins/
srun /export/users/sofia/sire-jm.app/bin/somd-freenrg -C ../../input/sim_ab.cfg -l $lam -p CUDA
wait
cd ..

