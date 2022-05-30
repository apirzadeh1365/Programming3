#!/bin/bash
#SBATCH --time 5:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH--partition=assemblix
#SBATCH --job-name=Azadeh
mkdir -p output
export BLASTDB=/local-fs/datasets/
export blastoutput=blastoutput.txt
export time=/usr/bin/time
export timings=output/timings.txt
for n in {1..16} ; do $time -a -o $timings -f %e blastp -query MCRA.faa -db ${BLASTDB}refseq_protein/refseq_protein -num_threads $n -outfmt 6 >> $blastoutput ; done
python assignment3.py