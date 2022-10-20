from Bio import SeqIO
import sys
dir = sys.argv[1]
list_of_kemer = []
with open("output/"+dir+"/contigs.fa") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seque = record.seq
        list_of_kemer.append([len(seque)])
list_of_kemer.sort(reverse=True, key=lambda x: x[0])
list_of_number = [item[0] for item in list_of_kemer]
avg = (sum(list_of_number))/2
sum_of_number = 0

for kmr in list_of_kemer:
    sum_of_number = sum_of_number + kmr[0]
    if sum_of_number >= avg:
        n50 = kmr[0]
        break
print(dir,n50)
