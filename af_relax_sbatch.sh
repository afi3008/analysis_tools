#!/bin/bash
#SBATCH -A NAISS2025-22-1521 -p alvis
#SBATCH --gpus-per-node=A40:1
#SBATCH -t 1-00:00:00
#SBATCH 

module load PyTorch-bundle/2.1.2-foss-2023a-CUDA-12.1.1
source /mimer/NOBACKUP/groups/naiss2025-23-641/afi/my_venv/bin/activate
nvidia-cuda-mps-control -d 
#export WANDB_API_KEY=
#wandb login
python aa_from_ca_mapping.py --mapping_smaller mapping_new.dat --mapping_bigger mapping_10.dat --output_prefix testing
#python af_relax_test.py --output_prefix relaxed
#cryosphere_train --experiment_yaml parameters.yaml 
#cryosphere_analyze --experiment_yaml parameters_test_aa.yaml --model cryoSPHERE/ckpt39.pt --segmenter cryoSPHERE/seg39.pt --output_path output_test_aa --generate_structures --all_atom
