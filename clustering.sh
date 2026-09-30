#!/bin/bash
#SBATCH -A NAISS2025-22-1521 -p alvis
#SBATCH --gpus-per-node=A100:1
#SBATCH -t 1-00:00:00
#SBATCH 

module load PyTorch-bundle/2.1.2-foss-2023a-CUDA-12.1.1
source /mimer/NOBACKUP/groups/naiss2025-23-641/afi/my_venv/bin/activate
#nvidia-cuda-mps-control -d 
#export WANDB_API_KEY=
#wandb login
#python combine.py --path_structures cryoSPHERE/analyze63/predicted_structures/ --output_path clustering/combine.pdb
#python clustering_test.py
python aa_from_ca.py --aa_structure ../fitted_structure_centered.pdb --output_prefix aa
#cryosphere_train --experiment_yaml parameters.yaml 
#cryosphere_analyze --experiment_yaml parameters_test_aa.yaml --model cryoSPHERE/ckpt39.pt --segmenter cryoSPHERE/seg39.pt --output_path output_test_aa --generate_structures --all_atom
