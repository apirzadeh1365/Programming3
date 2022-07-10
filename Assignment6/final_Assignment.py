"""
Author: Azadeh Pirzadeh
This code will predict the InterPro_accession of
long feature based on small features
input:  tsv file
output: The accuracy of model
"""
import dask.dataframe as dd
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from dask_ml import preprocessing, model_selection

PATH= "/data/dataprocessing/interproscan/all_bacilli.tsv"
# read the data
df = pd.read_csv(PATH, header=None, dtype={8: "object"}, sep="\t", usecols=[0, 1, 2, 6, 7, 11])
# change the pandas dataframe to dask dataframe
ddf = dd.from_pandas(df, npartitions=20, sort=False).reset_index()
# select and rename the useful columns
ddf = ddf[[0, 1, 2, 6, 7, 11]].rename(
    columns={0: "Protein_accession", 1: "MD5", 2: "Seq_len",
             6: "Start", 7: "Stop", 11: "InterPro_accession"})
# delete the "-"
ddf = ddf[ddf["InterPro_accession"] != "-"]
# make a column which contains the length of feature(Stop-start)
ddf["S-S"] = ddf["Stop"] - ddf["Start"]
# divided the length of feature on protein length
ddf["S-S/len"] = ddf["S-S"] / ddf["Seq_len"]


def class_p(data):
    """
    This function will find the long feature and small features
    input:  dask dataframe
    output: 0 for small and 1 for long
    """
    if data["S-S/len"] > 0.9:
        return 1
    else:
        return 0


# make a column which show the long and small feature(1:long and 0:small)
ddf["class_p"] = ddf.apply(class_p, meta=(None, 'int64'), axis=1)


def status(data):
    """
      This function will find proteins which have at
      least one small and one large features
      input:  dask dataframe
      output: True for useful protein and false foe unuseful
      """
    if 0 < np.mean(data['class_p']) < 1:
        return True
    else:
        return False


# find the list of the proteins which have at least one small
# and one large features and change it to list type
Protein = ddf.groupby("Protein_accession").\
    apply(lambda x: status(x), meta=('status', 'f8')).reset_index()
Protein_list = Protein[Protein['status'] == True]
Protein_list = Protein_list.compute()["Protein_accession"].to_list()
# keep the useful proteins from whole dataset
ddf = ddf[ddf["Protein_accession"].isin(Protein_list)]
# create 2 new dataset one with long features and another with small feature
ddf_long = ddf[ddf['class_p'] == 1]
ddf_small = ddf[ddf['class_p'] == 0]
# find the largest feature of each protein and save it in another dataframe
ddf_class = ddf_long.groupby(["Protein_accession", "InterPro_accession"])["S-S"]\
    .agg("max").reset_index()
# keep the important columns
ddf_class = ddf_class[["Protein_accession", "InterPro_accession"]]
# keep the important feature of small_feature dataset
ddf_small = ddf_small[["Protein_accession", "InterPro_accession"]]
# add one column to dataset with value of 1 to use in pivot tabel
ddf_small['count'] = 1
# categorize the InterPro_accession
ddf_small = ddf_small.categorize(columns=['InterPro_accession'])
# make a tabel of counts of small features of each protein
Matrix = dd.reshape.pivot_table(ddf_small, index="Protein_accession",
                                columns="InterPro_accession", values="count",
                                aggfunc='sum')
# merge the small features and large features to one dataframe
joined = Matrix.merge(ddf_class, how="left", on=["Protein_accession"])
# make a X and y to start the machine learning technic
X = joined.iloc[:, 1:-1]
y = joined[["InterPro_accession"]]
# categorize and use on hot encoding to predict just one protein
y = y.categorize()
y = preprocessing.OneHotEncoder().fit_transform(y)
# split the data to train and test
X_train, X_test, y_train, y_test = model_selection.train_test_split(X, y, test_size=0.3)
# make a random forest model
clf = RandomForestClassifier(n_estimators=200)
clf.fit(X_train, y_train)
y_pred = clf.predict(X_test)
# calculate and save the accuracy of model in a text file
Acuuracy=metrics.accuracy_score(y_test, y_pred)
with open('result.txt', 'w') as f:
    f.write("Accuracy : ")
    f.write(str(Acuuracy))