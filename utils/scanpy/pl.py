#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""utils.scanpy.pl
用于为adata绘图
"""


# In[ ]:


from utils.scanpy.sc import Path,np,pd,sc
from utils.plot import colormap


# In[ ]:


def spatial(adata, key, key_uns_spatial, key_img, scale_factor, ax,
                  color_map=None, spot_size=10, marker='.',
                  draw_img=True, draw_scatter=True,
                  kw_imshow={}
                  ):
    # 前处理， 判断与赋值
    assert key in adata.obs.columns
    assert key_uns_spatial in adata.uns['spatial']
    assert key_img in adata.uns['spatial'][key_uns_spatial]['images']
    assert draw_scatter or draw_img, '[Error] at least one of draw_img and draw_scatter must be True'
    if color_map is None:
        colormap.get(np.unique(adata.obs[key]))
    df_spatial = pd.DataFrame(adata.obsm['spatial'],
                              index=adata.obs.index,
                              columns='spatial1,spatial2'.split(','))\
        .mul(scale_factor)\
        .join(adata.obs.loc[:, [key]]).copy()
    df_spatial[key] = df_spatial[key].astype(str)
    dict_spatial = adata.uns['spatial'][key_uns_spatial]
    # scatter
    [ax.scatter('spatial1', 'spatial2', label=label, s=spot_size,
                marker=marker, c=color_map[label],
                data=df_spatial.query("{} == '{}'".format(key, label)))
        for label in color_map.keys()]
    img_edge = np.concatenate([ax.get_xlim(), ax.get_ylim()])
    ax.clear() if not draw_scatter else None
    # im show
    if draw_img:
        ax.imshow(dict_spatial['images'][key_img], **kw_imshow)
        ax.grid(False)
        ax.set_xlim(img_edge[0], img_edge[1])
        ax.set_ylim(img_edge[2], img_edge[3])

    if not ax.yaxis_inverted():  # 确保y轴是转置的，符合img的坐标习惯
        ax.invert_yaxis()
    ax.set_frame_on(False)
    ax.set_xticks([], []), ax.set_yticks([], [])

