from pyspark.sql import SQLContext
from pyspark import SparkContext
from pyspark.sql import functions as f
import pandas as pd
import numpy as np

sc=SparkContext('local[16]')
path= "/data/dataprocessing/interproscan/all_bacilli.tsv"
df = SQLContext(sc).read.csv(path, sep=r'\t', header=False, inferSchema= True)

df = df.withColumnRenamed('_c0', 'Protein_accession')
df = df.withColumnRenamed('_c1', 'MD5')
df = df.withColumnRenamed('_c2', 'Seq_len')
df = df.withColumnRenamed('_c3', 'Analysis')
df = df.withColumnRenamed('_c4', 'Signature_accession')
df = df.withColumnRenamed('_c5', 'Signature_description')
df = df.withColumnRenamed('_c6', 'Start')
df = df.withColumnRenamed('_c7', 'Stop')
df = df.withColumnRenamed('_c8', 'Score')
df = df.withColumnRenamed('_c9', 'Status')
df = df.withColumnRenamed('_c10', 'Date')
df = df.withColumnRenamed('_c11', 'InterPro_accession')
df = df.withColumnRenamed('_c12', 'InterPro_discription')
df = df.withColumnRenamed('_c13', 'GO_annotations')
df = df.withColumnRenamed('_c14', 'Pathways')

A1= df.filter(df.InterPro_accession!="-").select(f.countDistinct("InterPro_accession"))
E1=A1._sc._jvm.PythonSQLUtils.explainString(A1._jdf.queryExecution(),"simple")
A1=A1.collect()[0][0]

A2=df.filter(df.InterPro_accession!="-").groupBy('Protein_accession').agg(f.count('InterPro_accession')).agg(f.mean('count(InterPro_accession)'))
E2=A2._sc._jvm.PythonSQLUtils.explainString(A2._jdf.queryExecution(),"simple")
A2=A2.collect()[0][0]

A3=df.filter(df.GO_annotations!="-").withColumn('word', f.explode(f.split(f.col('GO_annotations'), '\|'))) \
  .groupBy('word').count().sort('count', ascending=False)
E3=A3._sc._jvm.PythonSQLUtils.explainString(A3._jdf.queryExecution(),"simple")
A3=A3.collect()[0][0]

A4=df.select(df.Stop-df.Start).agg(f.mean('(Stop - Start)'))
E4=A4._sc._jvm.PythonSQLUtils.explainString(A4._jdf.queryExecution(),"simple")
A4=A4.collect()[0][0]

A5=A5=df.filter(df.InterPro_accession!="-").groupBy("InterPro_accession").count().sort(f.desc("count"))
E5=A5._sc._jvm.PythonSQLUtils.explainString(A5._jdf.queryExecution(),"simple")
A5=A5.select('InterPro_accession').head(10)
A5=[data[0] for data in A5]

df=df.withColumn('S_S',(df.Stop-df.Start)/(df.Seq_len))
A6=df.filter(df.InterPro_accession!="-").select('InterPro_accession').filter((df.S_S)>0.9).sort(f.desc(df.Stop-df.Start))
E6=A6._sc._jvm.PythonSQLUtils.explainString(A6._jdf.queryExecution(),"simple") 
A6=A6.head(10)
A6=[data[0] for data in A6]

A7=df.filter(df.InterPro_discription!="-").withColumn('word', f.explode(f.split(f.col('InterPro_discription'), '\s|,'))) \
  .groupBy('word').count().sort('count', ascending=False)
E7=A7._sc._jvm.PythonSQLUtils.explainString(A7._jdf.queryExecution(),"simple") 
A7=A7.select('word').where('word != ""').head(10)
A7=[data[0] for data in A7]

A8=df.filter(df.InterPro_discription!="-").withColumn('word', f.explode(f.split(f.col('InterPro_discription'), '\s|,'))) \
  .groupBy('word').count().sort('count', descending=False)
E8=A8._sc._jvm.PythonSQLUtils.explainString(A8._jdf.queryExecution(),"simple") 
A8=A8.select('word').where('word != ""').head(10)
A8=[data[0] for data in A8]

A9=df.filter(df.InterPro_discription!="-").filter((df.S_S)>0.9).sort(f.desc(df.Stop-df.Start)).withColumn('word', f.explode(f.split(f.col('InterPro_discription'), '\s|,'))) \
  .groupBy('word').count().sort('count', ascending=False)
E9=A9._sc._jvm.PythonSQLUtils.explainString(A9._jdf.queryExecution(),"simple") 
A9=A9.select('word').where('word != ""').head(10)
A9=[data[0] for data in A9]

A10=df.filter(df.InterPro_accession !="-").groupBy('Protein_accession','Seq_len').count()
E10 = A10._sc._jvm.PythonSQLUtils.explainString(A10._jdf.queryExecution(),"simple") 
A10=A10.corr('Seq_len','count')**2

Az_data=pd.DataFrame(columns=['Question_number','Answer','Explain'])
Az_data.Question_number=np.arange(1,11)
Az_data.Answer=[A1,A2,A3,A4,A5,A6,A7,A8,A9,A10]
Az_data.Explain=[E1,E2,E3,E4,E5,E6,E7,E8,E9,E10]
Az_data.to_csv('output/assignment5.csv', index=False,header=True)
