"""
assignment 01
Azadeh Pirzadeh
"""
from multiprocessing import Pool, cpu_count
import argparse as ap
from Bio import Entrez

Entrez.email = "az.pirzadeh@gmail.com"
API_KEY = "b7989dc34851872fc7c8fffe0ba425979708"


def search(pm_id):
    results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc",
                                       LinkName="pubmed_pmc_refs", id=pm_id, api_key=API_KEY))
    references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
    return references[:10]


def write(references):
    handle = Entrez.efetch(db="pmc", id=references, rettype="XML", retmode="text", api_key=API_KEY)
    with open(f'output/{references}.xml', 'wb') as file:
        file.write(handle.read())


if __name__ == '__main__':
    argparse = ap.ArgumentParser(
        description="Script that downloads (default) 10 articles"
                    " referenced by the given PubMed ID concurrently.")
    argparse.add_argument("pubmed_id", action="store", type=str, nargs=1,
                          help="Pubmed ID of the article to harvest for references to download.")
    args = argparse.parse_args()
    ref = search(args.pubmed_id)
    c = cpu_count()
    with Pool(c) as p:
        p.map(write, ref)
