from Bio import SeqIO
import sys

dir = sys.argv[1]
list_of_kemer = []
with open("output/"+dir+"/contigs.fa") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        sequence = record.seq
        list_of_kemer.append([len(sequence)])
list_of_kemer.sort(reverse=True, key=lambda x:x[0])
list_of_number = [item[0] for item in list_of_kemer]
ave = (sum(list_of_number))/2

sum_number = 0
for kmer in list_of_kemer:
    sum_number = sum_number + kmer[0]
    if sum_number >= ave:
        n50 = kmer[0]
        break
print(dir,n50)
