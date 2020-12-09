import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import os

path = os.path.join(
    os.sep, 
    "home", 
    "andrew", 
    "Documents", 
    "MSc", 
    "wsmmGroup", 
    "Cardiomyopathies - including childhood onset", 
    "moduland", 
    )

bridgnessCentralityFile = os.path.join(path, "bridgness_centrality.csv")
df = pd.read_csv(bridgnessCentralityFile) 
df['moduland_betweenness'] = df['moduland_betweenness'].apply(lambda x: float(x))
df['moduland_centrality'] = df['moduland_centrality'].apply(lambda x: float(x))
# print(df)
# sorting modules into sub df's
modules = {
    "MLCYD":df[(df['module'] == "MLCYD")],
    "NTRK1":df[(df['module'] == "NTRK1")],
    "HCN1":df[(df['module'] == "HCN1")],
}

genes = df.name

panelAppGreen = df[(df['panelAppGreen'] == 1)]

# Matplotlib
x = df.moduland_centrality
y = df.moduland_betweenness

xMLCYD = modules["MLCYD"].moduland_centrality
yMLCYD = modules["MLCYD"].moduland_betweenness

xNTRK1 = modules["NTRK1"].moduland_centrality
yNTRK1 = modules["NTRK1"].moduland_betweenness

xHCN1 = modules["HCN1"].moduland_centrality
yHCN1 = modules["HCN1"].moduland_betweenness

xPanelAppGreen = panelAppGreen.moduland_centrality
yPanelAppGreen = panelAppGreen.moduland_betweenness


def regression(x,y):
    a, b = np.polyfit(np.array(x), np.array(y), deg=1)
    f = lambda x: a*np.array(x) + b
    return f

fMLCYD = regression(xMLCYD, yMLCYD)
fNTRK1 = regression(xNTRK1, yNTRK1)
fHCN1 = regression(xHCN1, yHCN1)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlabel("Centrality")
ax1.set_ylabel("Betweeness")
# ax1.scatter(x, y, c='k', marker="s", label='All entries')
ax1.scatter(xMLCYD, yMLCYD, c='r', marker="s", label='MLCYD')
ax1.scatter(xNTRK1, yNTRK1, c='gold', marker="s", label='NTRK1')
ax1.scatter(xHCN1, yHCN1, c='g', marker="s", label='HCN1')
ax1.scatter(xPanelAppGreen, yPanelAppGreen, c='k', marker="x", label='panelAppGreen')

ax1.plot(xMLCYD,fMLCYD(xMLCYD),lw=1.5, c="r")
ax1.plot(xNTRK1,fNTRK1(xNTRK1),lw=1.5, c="gold")
ax1.plot(xHCN1,fHCN1(xHCN1),lw=1.5, c="g")

plt.legend(loc='upper left')

# for i, label in enumerate(genes):
#     ax1.annotate(label, (x[i], y[i]))

plt.show()