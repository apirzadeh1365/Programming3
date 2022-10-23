#!/bin/bash
#SBATCH --job-name=Azadeh
#SBATCH --time 48:00:00
#SBATCH --nodes=10
#SBATCH --mem=1000
#SBATCH --cpus-per-task=16
#SBATCH --partition=assemblix
az_folder=output
mkdir -p ${az_folder}
Data1=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R1_001_BC24EVACXX.filt.fastq
Data2=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R2_001_BC24EVACXX.filt.fastq
output=${az_folder}/
seq 25 2 31 | parallel -j16 "velveth ${output}{} {} -longPaired -fastq -separate ${Data1} ${Data2}"
output=${az_folder}/
seq 25 2 31 | parallel -j16 "velvetg ${output}{}"
seq 25 2 31 | parallel python3 assignment4.py {} >> output/kmer_result.csv
python3 assignment5.py