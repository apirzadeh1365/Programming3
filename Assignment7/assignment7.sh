#!/bin/bash
#SBATCH --job-name=Azadeh
#SBATCH --time 48:00:00
#SBATCH --nodes=10
#SBATCH --mem=1000
#SBATCH --cpus-per-task=16
#SBATCH --partition=assemblix
mkdir -p output
path="/data/datasets/NCBI/genbank/Bacteria"
az_folder=output
find ${path} -type f -name "*.gbk" | parallel python3 assignment7.py {} >> ./all_composition.txt
