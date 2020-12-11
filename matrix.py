import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import os
import random
import seaborn as sns

# path to file
path = os.path.join(
    os.sep, 
    "media", 
    "sf_W_DRIVE", 
    "OneDrive - UOM", 
    "MSc", 
    "Year 3", 
    "wsmmGroup", 
    "Cardiomyopathies - including childhood onset", 
    )
matrixFile = os.path.join(path, "All_Genes_Matrix.csv")

df = pd.read_csv(matrixFile, index_col=0) 
dfTranspose = df.transpose()
# print(dfTranspose)

print(" ")

print(df.dot(dfTranspose))

geneMatrix = df.dot(dfTranspose)
# coreGeneMatrix.to_csv(path+os.sep+"coreGeneMatrix.csv")


hpoMatrix = dfTranspose.dot(df)
# coreHpoMatrix.to_csv(path+os.sep+"coreHpoMatrix.csv")

print(" ")

print(dfTranspose.dot(df))
# result = np.dot(df, dfTranspose)
# print(result)
# result= np.dot(df, dfTranspose)

# print(dfTranspose.dtypes)


# print(df.transpose())
df = hpoMatrix

plt.pcolor(df)
plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)
plt.show()
