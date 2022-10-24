"""
assignment4
Azadeh Pirzadeh
"""
from Bio import SeqIO
import sys
import numpy as np

dir = sys.argv[1]
list_of_kemer = []
with open("output/"+dir+"/contigs.fa") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        sequence = record.seq
        list_of_kemer.append(len(sequence))   
list_of_kemer.sort(reverse=True)
threshold = (np.sum(list_of_kemer))/2
sum_number = 0
for kmr in list_of_kemer:
    while sum_number <= threshold:
        sum_number = sum_number + kmr
        n50 = kmr
        break
print(dir,n50)
