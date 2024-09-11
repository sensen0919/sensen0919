import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import scrublet as scr
import scvi
import scib
import pymde
from scipy.stats import median_abs_deviation
import scanpy as sc
import omicverse as ov
ov.plot_set()
adata_T=sc.read("/Users/wenshidai/Desktop/python/adata_T.h5ad")
ov.pp.scale(adata_T)
ov.pp.pca(adata_T,layer='scaled',n_pcs=50)

adata_T.obsm["X_mde_pca"] = ov.utils.mde(adata_T.obsm["scaled|original|X_pca"])
adata_T=ov.single.batch_correction(adata_T,batch_key='sample',
                           methods='scVI',n_layers=2, n_latent=30, gene_likelihood="nb")
adata_T
adata_T.obsm["X_mde_scVI"] = ov.utils.mde(adata_T.obsm["X_scVI"])
sc.pp.neighbors(adata_T, n_neighbors=15, n_pcs=50,
               use_rep='scaled|original|X_pca')
sc.tl.leiden(
    adata_T,
    resolution=1.5,
    random_state=0,
    flavor="igraph",
    n_iterations=2,
    directed=False,
)
sc.tl.dendrogram(adata_T,'leiden',use_rep='scaled|original|X_pca')
sc.tl.rank_genes_groups(adata_T, 'leiden', use_rep='scaled|original|X_pca',
                        method='wilcoxon',use_raw=False,)
adata_T.obsm["X_mde"] = ov.utils.mde(adata_T.obsm["scaled|original|X_pca"])
mark=["KLRB1","ZBTB16","SLC4A10"]
mark1=["TRDC"]
mark2=["MKI67"]
mark3=["CD4","CD8A"]
mark4=["BCL6","CXCR5","PDCD1","CD4"]
mark5=["IL17A","IL17F"] #IL17
mark6=["IL2RA","FOXP3","CTLA4","TNFRSF4"]
sc.pl.dotplot(adata_T,mark3,groupby="leiden",cmap="RdBu_r",standard_scale="var")
Celltype={
    "0":"CD4",
    ...
    }
adata_T.obs["Celltype"]=adata_T.obs["leiden"].map(Celltype)
#adata_CD8&adata_CD4 same
adata_CD4=adata_CD4[(adata_CD4.obs['leiden']!='10')&(adata_CD4.obs['leiden']!='21')].copy()
adata_other2=adata_CD8[(adata_CD8.obs['leiden']=='24')|(adata_CD8.obs['leiden']=='25')].copy()
a=adata_CD8.obs['Celltype'].copy()
a=pd.DataFrame(a)
b=adata_CD4.obs["Celltype"].copy()
b=pd.DataFrame(b)
c=adata_DP.obs["Celltype"].copy()
c=pd.DataFrame(c)
d=adata_pro.obs["Celltype"].copy()
d=pd.DataFrame(d)
e=adata_Tgd.obs["Celltype"].copy()
e=pd.DataFrame(e)
f=adata_MAIT.obs["Celltype"].copy()
f=pd.DataFrame(f)
dfs=[a,b,c,d,e,f]
h=pd.concat(dfs)
