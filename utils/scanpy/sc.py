#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""sc scanpy


"""


# In[3]:


from utils.general import Path,np,pd,plt,mpl,sns,json
import scanpy as sc
import scipy
from collections.abc import Iterable

from PIL import Image
# 有的组织图像太大，取消Image的图像大小限制
Image.MAX_IMAGE_PIXELS = None


# # I/O

# In[ ]:


def load_adata(p_dir,prefix=''):
    def _load_json(p):
        return json.loads(p.read_text())

    def load_h5ad_from_mtx(p_dir,prefix=''):
        p_dir = Path(p_dir)
        assert p_dir.joinpath('{}matrix.mtx'.format(prefix)).exists(
        ), '[not exists] {}matrix.mtx\nin {}'.format(prefix,p_dir)
        
        adata = sc.read_10x_mtx(p_dir,prefix=prefix)
        # obs.csv
        if p_dir.joinpath('{}obs.csv'.format(prefix)).exists():
            adata.obs = pd.read_csv(p_dir.joinpath('{}obs.csv'.format(prefix)), index_col=0)
        else:
            print('[not exists]{}obs.csv\nin {}'.format(prefix,p_dir))
        return adata
    
    p_dir = Path(p_dir)
    adata = None
    if p_dir.match("*.h5ad"):
        adata = sc.read_h5ad(p_dir)
    elif p_dir.is_dir() and p_dir.joinpath('{}matrix.mtx'.format(prefix)).exists():
        adata = load_h5ad_from_mtx(p_dir,prefix)
    else:
        raise Exception("[can not load adata] {}".format(p_dir))

    # [load] spatial info: adata.obsm and adata.uns['spatial']
    ## adata.obsm
    if p_dir.joinpath('{}obsm'.format(prefix)).exists():
        for p_obsm in p_dir.joinpath('{}obsm'.format(prefix)).iterdir():
            if not p_obsm.match('*csv'):
                continue
            adata.obsm[p_obsm.stem] =pd.read_csv(p_obsm).to_numpy()
    ## adata.uns['spatial']
    if p_dir.joinpath('{}uns/spatial'.format(prefix)).exists():
        adata.uns['spatial'] = {}
        for p_uns_spatial in p_dir.joinpath('{}uns/spatial'.format(prefix)
                                           ).iterdir():
            
            dict_spatial = {'images':{}}
            for img in p_uns_spatial.joinpath('images').iterdir():
                dict_spatial['images'][img.stem]  = plt.imread(img)
                
            for p in p_uns_spatial.iterdir():
                if p.match('*.json'):
                    dict_spatial[p.stem] = _load_json(p)

            adata.uns['spatial'][p_uns_spatial.stem] = dict_spatial

    return adata

def load_spatial_images(adata,key_uns_spatial=None,
    key_images_path='images_path',update_images=False):
    def load_spatial_images_one(
            adata,
            key_uns_spatial,
            key_images_path='images_path',
            update_images=False):
        assert key_uns_spatial in adata.uns['spatial'].keys()
        dict_spatial = adata.uns['spatial'][key_uns_spatial]

        assert key_images_path in dict_spatial.keys()
        assert dict_spatial[key_images_path], "[Error] no img path"

        if 'images' not in dict_spatial.keys():
            dict_spatial['images'] = {}

        for k_img_path, img_path in dict_spatial[key_images_path].items():
            if update_images or (
                    k_img_path not in dict_spatial['images'].keys()):
                dict_spatial['images'][k_img_path] = mpl.image.imread(
                    img_path)

        adata.uns['spatial'][key_uns_spatial] = dict_spatial
        return adata

    assert 'spatial' in adata.uns.keys()
    if key_uns_spatial is None:
        key_uns_spatial = adata.uns['spatial'].keys()
    elif isinstance(key_uns_spatial, str):
        key_uns_spatial = [key_uns_spatial]

    for k in key_uns_spatial:
        adata = load_spatial_images_one(
            adata, k, key_images_path, update_images)
    return adata

def save_as_mtx(adata, p_dir, layer='counts', prefix='', as_int=True,
               key_images_path='images_path'):
    """
    将adata对象保存为matrix.mtx,barcodes.tsv,genes.tsv
    尝试保存obs.csv
    尝试保存空间信息,adata.obsm and adata.uns['spatial']
    p_dir ： 输出路径
    as_int : 是否将矩阵转换int类型
        default True
"""
    def _save_img(arr,p):
        # 存为npz和npy都极其巨大
        # 图片压缩技术是真的强
        im = Image.fromarray(arr)
        im.save(p)
    
    def _save_json(data,p):
        p.write_text(json.dumps(data))
    
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
        p_dir.joinpath("{}genes.tsv".format(prefix)),
        header=False,
        index=False,
        sep="\t",
    )
    # [out] barcodes.tsv obs.csv
    adata.obs.loc[:, []].to_csv(
        p_dir.joinpath("{}barcodes.tsv".format(prefix)),
        header=False,
        index=True,
        sep="\t",
    )

    if len(adata.obs_keys()) > 0:
        adata.obs.to_csv(
            p_dir.joinpath("{}obs.csv".format(prefix)), index=True
        )

    # [out] matrix.mtx
    data = adata.layers[layer] if layer in adata.layers.keys() else adata.X
    data = scipy.sparse.csr_matrix(data)
    if as_int:
        data = data.astype(int)
    nonzero_index = [i[:10] for i in data.nonzero()]
    print(
        "frist 10 matrix nonzero elements:\n",
        data[nonzero_index[0], nonzero_index[1]],
    )
    scipy.io.mmwrite(
        p_dir.joinpath("{}matrix.mtx".format(prefix)), data.getH()
    )

    # [save] spatial info: adata.obsm and adata.uns['spatial']
    if 'spatial' in adata.uns.keys() and 'spatial' in adata.obsm.keys():
        print('[save] spatial info')
        # [save] adata.obsm
        p_dir.joinpath(
            '{}obsm'.format(prefix)).mkdir(
            parents=True,
            exist_ok=True)
        for k_obsm in adata.obsm.keys():
            pd.DataFrame(
                adata.obsm[k_obsm],
                index=adata.obs.index,
                columns=[
                    '{}{}'.format(
                        k_obsm,
                        _+
                        1) for _ in range(
                        adata.obsm[k_obsm].shape[1])]) .to_csv(
                    p_dir.joinpath(
                        '{}obsm/{}.csv'.format(
                            prefix,
                            k_obsm)),
                index=False)

        # [save] adata.uns['spatial']
        for k_spatial, dict_spatial in adata.uns['spatial'].items():
            p_dir.joinpath('{}uns/spatial/{}/images'.format(prefix,
                           k_spatial)) .mkdir(parents=True, exist_ok=True)
            # img 仅保存dict_spatial[key_images_path]中没有路径的img
            
            for k_img, img in dict_spatial.setdefault('images',{}).items():
                if k_img in  dict_spatial.setdefault(key_images_path,{}).keys():
                    pass
                else:
                    _save_img(img,p_dir.joinpath(
                        '{}uns/spatial/{}/images/{}.jpg'\
                        .format(prefix, k_spatial, k_img)))
            
            if dict_spatial.setdefault(key_images_path,None):
                _save_json({_k : str(Path(_v).absolute())
                for _k,_v in dict_spatial[key_images_path].items()},
                           p_dir.joinpath('{}uns/spatial/{}/{}.json'\
                                          .format(prefix, k_spatial,key_images_path)))
            # scalefactors
            if dict_spatial.setdefault('scalefactors',None):
                _save_json(dict_spatial['scalefactors'],
                           p_dir.joinpath('{}uns/spatial/{}/scalefactors.json'\
                                          .format(prefix, k_spatial)))
            # metadata
            if dict_spatial.setdefault('metadata',None):
                _save_json(dict_spatial['metadata'],
                           p_dir.joinpath('{}uns/spatial/{}/metadata.json'\
                                          .format(prefix, k_spatial)))

    print("[out] {}".format(p_dir))


# In[ ]:


def qc(adata):
    adata.var["mt"] = adata.var_names.str.match('^MT-',case=False)
    adata.var["ribo"] = adata.var_names.str.match('^RP[SL]',case=False)
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]",case=False)
    sc.pp.calculate_qc_metrics(adata, qc_vars='mt,ribo,hb'.split(','),percent_top=[],
                               log1p=False,inplace=True)


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

def subset_adata(adata, *args):
    def _process_values(values):
        if isinstance(values, Iterable):
            if isinstance(values, str):
                values = [values]
        else:
            values = [values]
        return values
    assert len(
        args) % 2 == 0, '[Error][{}] length of args must be 2*n'.format(len(args))

    for key, values in zip(args[0::2], args[1::2]):
        values = _process_values(values)
        adata = adata[adata.obs[key].isin(values), :]
    return adata

