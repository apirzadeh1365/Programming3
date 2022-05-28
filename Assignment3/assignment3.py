import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
df=pd.read_csv('output/timings.txt',header=None)
x=np.arange(1,17)
y=df[0][-16:].astype('float')
print(y.dtype)
plt.plot(x,y)
plt.xlabel('cpu')
plt.ylabel('Time')
plt.title("timing plot")
plt.savefig("output/timings.png")
print("Hooooray, It is done. You can see my plot")