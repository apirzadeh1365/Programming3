"""
Second part of assignment4
Azadeh Pirzadeh
"""
import pandas as pd
import numpy as np
import shutil as sh


df=pd.read_csv('output/kmer_result.csv',sep=' ',names=['kmr','n50'])
i=np.max(df['n50'])
number_i=df[df['n50']==i]
best=number_i.iloc[0,0]
source= r'output/{}/contigs.fa'.format(best)
destination= r"output/"
sh.copy(source,destination)
del_list = list(df['kmr'])
for i in del_list:
    sh.rmtree("output/"+str(i))