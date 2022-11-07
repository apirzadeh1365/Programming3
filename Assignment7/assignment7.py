"""
assignment6
Azadeh Pirzadeh
"""
import sys
from Bio import SeqIO
import pathlib
# data="/data/datasets/NCBI/genbank/Bacteria/Nitrosococcus_watsoni_C_113_uid40331/CP002086.gbk"
# data1 = "/data/datasets/NCBI/genbank/Bacteria/Nitrosococcus_watsoni_C_113_uid40331/*.gbk"
# data_path="/data/datasets/NCBI/genbank/Bacteria/*/*.gbk"
output="/students/2021-2022/master/Azadeh_DSLS/output"

def count_amino(genbank_path):
    amino_acid = ['C', 'D', 'S', 'Q', 'K', 'P', 'T', 'F', 'A', 'X', 'G', 'I', 'E', 'L', 'H', 'R', 'W', 'M', 'N', 'Y', 'V']
    counts_all = {}
    counts_first = {}
    counts_last = {}
    for amino in amino_acid: counts_all[amino] = 0
    for amino in amino_acid: counts_first[amino] = 0
    for amino in amino_acid: counts_last[amino] = 0

    total = 0
    total_fifth = 0
    text = ""
    try:
        for sequence in SeqIO.parse(genbank_path, "genbank"):
            for feature in sequence.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:

                    proteins=feature.qualifiers['translation']
                    for protein in proteins:
                        for a in amino_acid:
                            counts_all[a] += protein.count(a)
                            counts_first[a] += protein[:5].count(a)
                            counts_last[a] += protein[-5:].count(a)
                            total_fifth += 5
                            total += len(protein)
    
        line_all = ""
        line_first = ""
        line_last = ""
        for amino in amino_acid:
            line_all += f'{amino} :{round(float(counts_all[amino] / total  if total > 0 else 0) * 100, 2)} '
            line_first += f'{amino} :{round(float(counts_first[amino] / total_fifth if total_fifth > 0 else 0) * 100,2)} '
            line_last += f'{amino} :{round(float(counts_last[amino] / total_fifth if total_fifth > 0 else 0) * 100 ,2)} '
        text=line_all+ "\n" +line_first + "\n"+ line_last+"\n"
    except:
        print(genbank_path)

    return text

def determine_protein(genbank_path):
    output_path=output +"/" +(pathlib.Path(genbank_path).name).replace(".gbk", "_composition.txt")
    txt=count_amino(genbank_path)
   
    with open(output_path,"w") as handle:
        handle.write(txt)


if len(sys.argv)>0:   
    genbank_file_path = sys.argv[1]
    print(genbank_file_path)
    determine_protein(genbank_file_path)
