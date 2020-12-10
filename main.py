import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import os
import random

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
    "moduland", 
    )
bridgnessCentralityFile = os.path.join(path, "bridgness_centrality.csv")

# set up pandas DF
df = pd.read_csv(bridgnessCentralityFile) 
df['moduland_betweenness'] = df['moduland_betweenness'].apply(lambda x: float(x))
df['moduland_centrality'] = df['moduland_centrality'].apply(lambda x: float(x))

# Pd series of gene names
genes = df.name

def regression(x,y):
    a, b = np.polyfit(np.array(x), np.array(y), deg=1)
    f = lambda x: a*np.array(x) + b
    return f

def plotModule(df, moduleName, plotAxes, colour, marker="s"):

    # Set variables
    x = df.moduland_centrality
    y = df.moduland_betweenness
    subDf = df[(df['module'] == moduleName)]
    xsubDf = subDf.moduland_centrality
    ysubDf = subDf.moduland_betweenness
    fsubDf = regression(xsubDf,ysubDf)

    # Plot chart
    plotAxes.scatter(xsubDf, ysubDf, c=colour, marker=marker, label=moduleName)
    # Plot regression line
    plotAxes.plot(xsubDf,fsubDf(xsubDf),lw=1.5, c=colour)

    # Annotate points
    lstGenes=subDf.name.tolist() 
    for i, gene in enumerate(genes):
        if gene in lstGenes:
            ax1.annotate(gene, (x[i], y[i]))

# Set-up plot
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlabel("Centrality")
ax1.set_ylabel("Betweeness")

# Plot moduland clusters (modules)
plotModule(df,"MLCYD", ax1, "royalblue")
plotModule(df,"NTRK1", ax1, "gold")
plotModule(df,"HCN1", ax1, "r")

# sub DF for panelApp Genes (green)
panelAppGreen = df[(df['panelAppGreen'] == 1)]
# plot panelApp Genes (Green)
x = df.moduland_centrality
y = df.moduland_betweenness
xPanelAppGreen = panelAppGreen.moduland_centrality
yPanelAppGreen = panelAppGreen.moduland_betweenness
ax1.scatter(xPanelAppGreen, yPanelAppGreen, c='limegreen', marker="x", label='panelAppGreen')
ax1.plot(xPanelAppGreen,regression(xPanelAppGreen,yPanelAppGreen)(xPanelAppGreen),lw=1.5, c="limegreen")
panelAppGreenGenes=panelAppGreen.name.tolist() 
for i, gene in enumerate(genes):
    if gene in panelAppGreenGenes:
        ax1.annotate(gene, (x[i], y[i]))

# Add legend
plt.legend(loc='upper left')

# plt.show()

##########################################

def extractHpoTerms(GeneSymbol):
    "Extract HPO terms from HPO API with a given gene symbol"
    import urllib.request
    import json

    entrezFromGene = "https://hpo.jax.org/api/hpo/search/?q={}" # gene symbol
    hpoFromEntrez ="https://hpo.jax.org/api/hpo/gene/{}" # entrez id

    data = json.load(urllib.request.urlopen(entrezFromGene.format(GeneSymbol)))
    entrezId = data["genes"][0]["entrezGeneId"] 
    data = json.load(urllib.request.urlopen(hpoFromEntrez.format(entrezId)))

    ontology = {}
    for entry in data["termAssoc"]:
        ontology[entry["name"]] = entry["ontologyId"]

    return ontology





class Module:

    def __init__(self, gene):
        self.moduleGene = gene
        self.geneObjects = ""

class Gene:

    def __init__(self, gene):
        self.gene = gene
        self.hpoTerms = ""
        self.hpoTermsDf = ""

    def setHpoTerms(self):
        try:
            self.hpoTerms = extractHpoTerms(self.gene)
        except:
            print("Unable to extract HPO terms for {}".format(self.gene))
            self.hpoTerms = {}

    def setHpoDf(self):

        data = {
            "phenotype": self.hpoTerms.keys(),
            "hpoId": self.hpoTerms.values(),
        }


        self.hpoTermsDf = pd.DataFrame.from_dict(data)


def setHpoClasses(df, moduleGene):

    module = Module(moduleGene)

    dfHPO = df[(df['module'] == moduleGene)]
    hpoNames = dfHPO.name.to_list()

    module.geneObjects = [Gene(i) for i in hpoNames]

    return module

module = setHpoClasses(df, "NTRK1")

for item in module.geneObjects:
    item.setHpoTerms()
    item.setHpoDf()
    print(item.gene)
    print(item.hpoTermsDf)
    print(" ")