#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""sc scanpy
用于操作pandas.DataFrame
> 
gourp_add
iterdir

"""


# In[3]:


from pathlib import Path
import numpy as np
import pandas as pd

import scanpy as sc

from scipy.io import mmwrite
from scipy.sparse import csr_matrix


# In[ ]:


def qc(adata):
    adata.var["mt"] = adata.var_names.str.match('^MT-',case=False)
    adata.var["ribo"] = adata.var_names.str.match('^RP[SL]',case=False)
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]",case=False)
    sc.pp.calculate_qc_metrics(adata, qc_vars='mt,ribo,hb'.split(','),percent_top=[],
                               log1p=False,inplace=True)


# In[4]:


def h5ad_to_mtx(adata, p_dir, prefixes="", as_int=True):
    """
    将adata对象保存为mtx
    p_dir ： 输出路径
    as_int : 是否将矩阵转换int类型
        default True
    """
    assert adata.obs.index.is_unique, '[Error] obs index is not unique'
    assert adata.var.index.is_unique, '[Error] var index is not unique'

    p_dir = Path(p_dir)
    p_dir.mkdir(parents=True, exist_ok=True)

    # [out] genes.tsv
    adata.var["gene_names"] = adata.var_names.to_numpy()
    if "gene_ids" not in adata.var_keys():
        adata.var["gene_ids"] = adata.var["gene_names"]
    df_genes = adata.var.loc[:, ["gene_ids", "gene_names"]]
    df_genes.to_csv(
        p_dir.joinpath("{}genes.tsv".format(prefixes)),
        header=False,
        index=False,
        sep="\t",
    )

    # [out] barcodes.tsv obs.csv
    adata.obs.loc[:, []].to_csv(
        p_dir.joinpath("{}barcodes.tsv".format(prefixes)),
        header=False,
        index=True,
        sep="\t",
    )

    if len(adata.obs_keys()) > 0:
        adata.obs.to_csv(
            p_dir.joinpath("{}obs.csv".format(prefixes)), index=True
        )

    # [out] matrix.mtx
    adata.X = csr_matrix(adata.X)
    if as_int:
        adata.X = adata.X.astype(int)
    nonzero_index = [i[:10] for i in adata.X.nonzero()]
    print(
        "frist 10 adata.X nonzero elements:\n",
        adata.X[nonzero_index[0], nonzero_index[1]],
    )
    mmwrite(
        p_dir.joinpath("{}matrix.mtx".format(prefixes)), adata.X.getH()
    )
    print("[out] {}".format(p_dir))


# In[ ]:


def standard_process(adata,kv_hvg = {'n_top_genes':2000,'batch_key':None},
    kv_neighbors = {'n_neighbors':15,'n_pcs':15},
    kv_umap = {},kv_leiden = {'resolution':.5}):
    """
单细胞标准流程

+ 将原始count矩阵存于adata.layers["counts"]
+ 标准化
    + sc.pp.normalize_total
    + sc.pp.log1p
+ hvgs
    + highly_variable_genes
+ 降维 PCA
    + sc.tl.pca
        + sc.pl.pca_variance_ratio
+ 降维 UMAP
    + sc.pp.neighbors
    + sc.tl.umap
    + sc.tl.leiden
        + sc.pl.umap
"""
    adata = adata.copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,**kv_hvg)
    sc.tl.pca(adata,n_comps=50)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
    sc.pp.neighbors(adata,**kv_neighbors)
    sc.tl.umap(adata,**kv_umap)
    sc.tl.leiden(adata,**kv_leiden)
    sc.pl.umap(adata, color=["leiden"])
    return adata

