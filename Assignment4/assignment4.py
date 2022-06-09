from Bio import SeqIO
import sys

dir = sys.argv[1]
kmer_list = []
with open("output/"+dir+"/contigs.fa") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        sequence = record.seq
        kmer_list.append([len(sequence)])
kmer_list.sort(reverse=True, key=lambda x: x[0])
num_list = [item[0] for item in kmer_list]
average = (sum(num_list))/2
sum_num = 0

for kmer in kmer_list:
    sum_num = sum_num + kmer[0]
    if sum_num >= average:
        n50 = kmer[0]
        break
print(dir,n50)
