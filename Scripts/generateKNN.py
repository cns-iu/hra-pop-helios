from pathlib import Path
from tqdm.auto import tqdm
import igraph as ig
from pynndescent import NNDescent
import pandas as pd
import numpy as np
import sklearn.decomposition
import xnetwork as xn
import scanpy as sc
import anndata
# Start by importing logger
import logging
logging.basicConfig(level=logging.INFO)



# Load the data
dataName = "GTEx_8_tissues_snRNAseq_atlas_071421.public_obs"
dataFolder = Path("Data")
dataPath = dataFolder / f"{dataName}.h5ad"
networksPath = dataFolder/"Networks"
networksPath.mkdir(exist_ok=True,parents=True)

# save log to file
logPath = dataFolder/f"{dataName}_preprocess.log"
fileHandler = logging.FileHandler(logPath)
logging.getLogger().addHandler(fileHandler)

logging.info(f"Processing {dataName}")
logging.info(f"Data path: {dataPath}")
logging.info(f"Networks path: {networksPath}")



enablePCA = True
useTruncatedSVD = False
useLogNormalization = True
kneighbors=10
pcaDimensions = 100
# sampleRate = 1.0

# log parameters

logging.info(f"Parameters:")
logging.info(f"PCA: {enablePCA}")
logging.info(f"Log normalization: {useLogNormalization}")
logging.info(f"kneighbors: {kneighbors}")
logging.info(f"PCA dimensions: {pcaDimensions}")
# logging.info(f"Sample rate: {sampleRate}")


adata = anndata.read_h5ad(dataPath)
list_cells = adata.obs.index.tolist()

networkProperties = [
#  'n_genes',
#  'fpr',
 'tissue',
#  'prep',
#  'individual',
#  'nGenes',
#  'nUMIs',
 'PercentMito',
 'PercentRibo',
 'Age_bin',
 'Sex',
#  'Sample ID',
#  'Participant ID',
#  'Sample ID short',
#  'RIN score from PAXgene tissue Aliquot',
#  'RIN score from Frozen tissue Aliquot',
 'Autolysis Score',
#  'Sample Ischemic Time (mins)',
#  'Tissue Site Detail',
#  'scrublet',
#  'scrublet_score',
#  'barcode',
#  'batch',
#  'n_counts',
#  'tissue-individual-prep',
 'Broad cell type',
 'Granular cell type',
#  'introns',
#  'junctions',
#  'exons',
#  'sense',
#  'antisense',
#  'intergenic',
#  'batch-barcode',
#  'exon_ratio',
#  'intron_ratio',
#  'junction_ratio',
 'log10_nUMIs',
 'leiden',
 'leiden_tissue',
#  'Tissue composition',
 'Cell types level 2',
 'Cell types level 3',
 'Broad cell type numbers',
 'Broad cell type (numbers)',
 'Tissue',
#  'channel'
 ]

# log used network properties
logging.info(f"Network properties:")
for prop in networkProperties:
    logging.info(prop)


logging.info(f"Dataset:")
logging.info(f"Number of cells: {adata.n_obs}")
logging.info(f"Number of genes: {adata.n_vars}")

# filter cells
logging.info(f"Filtering cells by highly variable genes")

adata.X = adata.layers["counts"]

sc.pp.highly_variable_genes(adata, flavor = "seurat_v3", n_top_genes=10000)
variableGeneIndices = adata.var.highly_variable
adata = adata[:, adata.var.highly_variable]



if(useLogNormalization):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
else:
    sc.pp.normalize_total(adata, target_sum=1e4)


if enablePCA:
    if useTruncatedSVD:
        pca = sklearn.decomposition.TruncatedSVD(n_components=pcaDimensions)
        # no need to make dense if truncatedSVD is used
        pca.fit(adata.X)
        networkMatrix = pca.transform(adata.X)
        # log importance % of variance explained
    else:
        pca = sklearn.decomposition.PCA(n_components=pcaDimensions)
        # convert to an array: this is necessary for the PCA

        pcaMatrix = adata.X.toarray()
        pcaMatrix = sklearn.preprocessing.StandardScaler().fit_transform(pcaMatrix)
        # make nan 0
        pcaMatrix[np.isnan(pcaMatrix)] = 0
        pca.fit(pcaMatrix)
        networkMatrix = pca.transform(pcaMatrix)
        # log importance % of variance explained
        logging.info("PCA explained variance ratio")
        # log cumulative importance
        logging.info(f"PCA cumulative importance: {np.cumsum(pca.explained_variance_ratio_)}")
else:
    networkMatrix = adata.X.copy()


variant = ""
if enablePCA:
    variant = "_PCA"
if useLogNormalization:
    variant += "_log"
# if(networksPath/f"knn_{kneighbors}_{variant}_pos.xnet").exists():
#     print(f"skipping knn_{kneighbors}_{variant}_pos.xnet")
# try 4 times with no exception
trials =4
while trials>0:
    try:
        knnData = NNDescent(networkMatrix,
                            metric="correlation",
                            n_neighbors=kneighbors,
                            verbose=True
        ).neighbor_graph
        break
    except Exception as e:
        print(e)
        trials -= 1
        if trials==0:
            break

if(trials==0):
    print(f"Failed to compute knn for {kneighbors} neighbors")

edges = []
weights = []
neighborIndices = []
for sourceIndex,theNeighbors in enumerate(tqdm(knnData[0])):
    for neighborIndex,neighbor in enumerate(theNeighbors):
        if(sourceIndex==neighbor):
            continue
        # correlation between source and neighbor
        correlationSimilarity = np.corrcoef(networkMatrix[sourceIndex,:],networkMatrix[neighbor,:])[0,1]
        edges.append((sourceIndex,neighbor))
        weights.append(correlationSimilarity)
        neighborIndices.append(neighborIndex)
    # edges.append((kneighbors[0],kneighbors[1]))
    # edges.append((kneighbors[0],kneighbors[2]))



g = ig.Graph(networkMatrix.shape[0],edges=edges,edge_attrs={"weight":weights, "neighborIndex":neighborIndices})
# remove any self loops
g = g.simplify(multiple=False,loops=True,combine_edges={"weight":"sum"})

g.vs["Label"] = list_cells
g.vs["Organ"] = [x.split("-")[-1] for x in list_cells]
for attributeName in networkProperties:
    if(attributeName in adata.obs.columns):
        g.vs[attributeName] = adata.obs[attributeName].tolist()
# g.vs["CellType"] = np.array(cellTypes)[sampleIndices]

logging.info(f"Graph properties:")
logging.info(f"Number of nodes: {g.vcount()}")
logging.info(f"Number of edges: {g.ecount()}")
xn.igraph2xnet(g,(networksPath/f"knn_{kneighbors}_{variant}.xnet"))
#remove links with 0 or negative weights
g.es.select(weight_lt=0).delete()
xn.igraph2xnet(g,(networksPath/f"knn_{kneighbors}_{variant}_pos.xnet"))

logging.info(f"Saved knn_{kneighbors}_{variant}.xnet and knn_{kneighbors}_{variant}_pos.xnet")
logging.info(f"Finished processing {dataName}")






