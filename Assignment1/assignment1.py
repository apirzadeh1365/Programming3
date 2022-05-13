from asyncore import read
from msilib import type_valid
import sys
from Bio import Entrez
from multiprocessing import Pool, cpu_count

def search(pmid):
   Entrez.email = "az.pirzadeh@gmail.com"
   results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs",id=pmid, 
                         api_key = "b7989dc34851872fc7c8fffe0ba425979708"))
   references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
   return references[:10]

def write(references):
   for ref in references:
    handle = Entrez.efetch(db="pmc", id=ref, rettype="XML", retmode="text", api_key = "b7989dc34851872fc7c8fffe0ba425979708")
    with open(f'E:\programming3\Assignment1\output\{ref}.xml', 'wb') as file:
         file.write(handle.read())


if __name__ == '__main__':
   
   pmid =str(sys.argv[1])
   refrence=search(pmid)
   c=cpu_count()
   with Pool(c) as p:
        p.map(write(refrence), refrence)


  