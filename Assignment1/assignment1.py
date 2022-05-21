from asyncore import read
from msilib import type_valid
import sys
from Bio import Entrez
import argparse as ap
from multiprocessing import Pool, cpu_count

Entrez.email = "az.pirzadeh@gmail.com"
def search(pmid):
   
   results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs",id=pmid, 
                         api_key = "b7989dc34851872fc7c8fffe0ba425979708"))
   references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
   return references[:10]

def write(ref):
   handle = Entrez.efetch(db="pmc", id=ref, rettype="XML", retmode="text", api_key = "b7989dc34851872fc7c8fffe0ba425979708")
   with open(f'E:\programming3\Assignment1\output\{ref}.xml', 'wb') as file:
         file.write(handle.read())


if __name__ == '__main__':

    argparser = ap.ArgumentParser(description="Script that downloads (default) 10 articles referenced by the given PubMed ID concurrently.")
   
    argparser.add_argument("pubmed_id", action="store", type=str, nargs=1, help="Pubmed ID of the article to harvest for references to download.")
    args = argparser.parse_args()
    refrence=search(args.pubmed_id)
    c=cpu_count()
    with Pool(c) as p:
        p.map(write, refrence)


  