#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""utils.scanpy.pl
用于为adata绘图
"""


# In[ ]:


from utils.scanpy.sc import Path, np, pd, plt, mpl, sns
from utils.scanpy.sc import sc
from utils.scanpy.sc import get_scalefactor,get_spot_size
import utils.plot as pl


# # umap

# In[ ]:


def umap(adata,key,ax,colormap=None,marker='.',size = 10,kw_scatter={}):
    assert key in adata.obs.columns
    if colormap is None:
            colormap = pl.colormap.get(adata.obs[key].sort_values().unique())
    df_umap = pd.DataFrame(adata.obsm['X_umap'],
                 index=adata.obs.index,
                 columns='UMAP1,UMAP2'.split(','))\
            .join(adata.obs.loc[:, [key]]).copy()
    df_umap[key] = df_umap[key].astype(str)
    # scatter
    [ax.scatter('UMAP1', 'UMAP2', label=label, s=size,
                marker=marker, c=colormap[label],
                data=df_umap.query("{} == '{}'".format(key, label)),
                **kw_scatter)
        for label in colormap.keys()]
    ax.set_axis_off()

def umap_gene(adata,key,ax,marker = '.',
    vmin=None, vmax=None,
    cmap='Reds',size=10,draw_cbar = True,
    kw_scatter={},
    kw_cbar = {'location':'right',
        'fraction':.025, # 在对应方向上的宽度与父轴宽度的比值
        'aspect':40,# cbar的长宽比
        'pad':.02,   # 与父轴的间距
        'format':'{x:.1f}'}):
    assert key in adata.obs.columns or key in adata.var_names
    
    df_umap = pd.DataFrame(adata.obsm['X_umap'],
                 index=adata.obs.index,
                 columns='UMAP1,UMAP2'.split(','))\
            .join(sc.get.obs_df(adata,[key])).copy()
    
    # scatter
    cbar = ax.scatter('UMAP1', 'UMAP2', s=size,c=df_umap[key],
                    marker=marker,vmin=vmin, vmax=vmax,cmap=cmap,
                    data=df_umap,**kw_scatter)
    ax.figure.colorbar(cbar,ax=ax,**kw_cbar) if draw_cbar else None
    ax.set_axis_off()
    
    return cbar


# # spatial

# In[ ]:


def spatial(adata,key,ax,key_uns_spatial='spatial',
        key_img='img',colormap=None,size=1,
        spot_size=None,scale_factor=None,
        marker='.',draw_img=True,draw_scatter=True,
        kw_scatter={},kw_imshow={},kw_get_scalefactor={}):
    """由空间信息和图像信息绘制图形, key对应的数据为离散变量

Parameters
----------
scale_factor : int,float (default : None)
    若为None, 尝试由util.scanpy.sc.default get_scalefactor()获取
    该函数默认值为1

size,spot_size : int,float (default : 1,None)
    控制scatter点的大小

    spot_size为None, 尝试由util.scanpy.sc.default get_spot_size()获取
    该函数默认值为1

    最后将size * scale_factor * spot_size * 0.5
    传入 ax.scatter
"""

    # 前处理， 判断与赋值
    assert key in adata.obs.columns
    assert key_uns_spatial in adata.uns['spatial']
    assert key_img in adata.uns['spatial'][key_uns_spatial]['images']
    assert draw_scatter or draw_img, '[Error] at least one of draw_img and draw_scatter must be True'
    if colormap is None:
        colormap = pl.colormap.get(adata.obs[key].sort_values().unique())
    if scale_factor is None:
        scale_factor = get_scalefactor(
            adata,
            key_uns_spatial=key_uns_spatial,
            key_img=key_img,
            **kw_get_scalefactor)
    if spot_size is None:
        spot_size = get_spot_size(adata, key_uns_spatial)
    circle_radius = size * scale_factor * spot_size * 0.5

    df_spatial = pd.DataFrame(adata.obsm['spatial'],
                              index=adata.obs.index,
                              columns='spatial1,spatial2'.split(','))\
        .mul(scale_factor)\
        .join(adata.obs.loc[:, [key]]).copy()
    df_spatial[key] = df_spatial[key].astype(str)
    dict_spatial = adata.uns['spatial'][key_uns_spatial]
    # scatter
    [ax.scatter('spatial1', 'spatial2', label=label, s=circle_radius,
                marker=marker, c=colormap[label],
                data=df_spatial.query("{} == '{}'".format(key, label)),
                **kw_scatter)
        for label in colormap.keys()]
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
    ax.set_axis_off()

def spatial_gene(
    adata, key, key_uns_spatial, key_img, scale_factor, ax,
    vmin=None, vmax=None,cmap='Reds', size=1,
    spot_size=None, marker='.',draw_img=True,
    draw_scatter=True, draw_cbar=True,
    kw_scatter={},kw_imshow={},
    kw_cbar={
        'location': 'right',
        'fraction': .025,  # 在对应方向上的宽度与父轴宽度的比值
        'aspect': 40,  # cbar的长宽比,fraction*aspect 即为cbar的长的占比
        'pad': .02,   # 与父轴的间距
        # 'format': '{x:.1f}'
    }
):
    """由空间信息和图像信息绘制图形, key对应的数据为连续性变量

Parameters
----------
scale_factor : int,float (default : None)
    若为None, 尝试由util.scanpy.sc.default get_scalefactor()获取
    该函数默认值为1

size,spot_size : int,float (default : 1,None)
    控制scatter点的大小

    spot_size为None, 尝试由util.scanpy.sc.default get_spot_size()获取
    该函数默认值为1

    最后将size * scale_factor * spot_size * 0.5
    传入 ax.scatter
"""
    # 前处理， 判断与赋值
    assert key in adata.obs.columns or key in adata.var_names
    assert key_uns_spatial in adata.uns['spatial']
    assert key_img in adata.uns['spatial'][key_uns_spatial]['images']
    assert draw_scatter or draw_img, '[Error] at least one of draw_img and draw_scatter must be True'
    if scale_factor is None:
        scale_factor = get_scalefactor(adata,key_uns_spatial=key_uns_spatial,
                                     key_img = key_img,**kw_get_scalefactor)
    if spot_size is None:
        spot_size = get_spot_size(adata, key_uns_spatial)
    circle_radius = size * scale_factor * spot_size * 0.5
    
    df_spatial = pd.DataFrame(adata.obsm['spatial'],
                              index=adata.obs.index,
                              columns='spatial1,spatial2'.split(','))\
        .mul(scale_factor)\
        .join(sc.get.obs_df(adata, [key]))

    dict_spatial = adata.uns['spatial'][key_uns_spatial]

    # scatter
    cbar = ax.scatter('spatial1', 'spatial2',
                      s=circle_radius, marker=marker, c=df_spatial[key],
                      cmap=cmap, vmax=vmax, vmin=vmin,
                      data=df_spatial, **kw_scatter)
    ax.figure.colorbar(cbar, ax=ax, **kw_cbar) if draw_cbar else None

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
    ax.set_axis_off()

    return cbar

def spatial_3d(adata,key,ax,colormap=None,scale_factor = None,
    marker = '.',size = 10,kw_scatter = {},height = 5,
    query_3d_line = '',kw_line = {'linewidth': .5, 'color': 'grey'},
    kw_view_init={'elev': 35, 'azim': -90, 'roll': 0}):
    """
Parameters
----------
scale_factor : int,float (default : None)
    若为None, 尝试由util.scanpy.sc.default get_scalefactor()获取
    该函数默认值为1

size,spot_size : int,float (default : 1,None)
    控制scatter点的大小

    spot_size为None, 尝试由util.scanpy.sc.default get_spot_size()获取
    该函数默认值为1

    最后将size * scale_factor * spot_size * 0.5
    传入 ax.scatter
"""
    assert key in adata.obs.columns
    if colormap is None:
        colormap = pl.colormap.get(adata.obs[key].sort_values().unique())
    if scale_factor is None:
        scale_factor = get_scalefactor(adata,key_uns_spatial=key_uns_spatial,
                                     key_img = key_img,**kw_get_scalefactor)
    if spot_size is None:
        spot_size = get_spot_size(adata, key_uns_spatial)
    circle_radius = size * scale_factor * spot_size * 0.5
    
    df_plot = pd.DataFrame(adata.obsm['spatial'],
                           index=adata.obs.index, columns='spatial1,spatial2'.split(','))\
        .mul(scale_factor).join(pd.DataFrame(adata.obsm['X_umap'],
                                             index=adata.obs.index, columns='UMAP1,UMAP2'.split(',')))\
        .join(adata.obs.loc[:, [key]])
    df_plot[key] = df_plot[key].astype(str)
    
    from utils.arr import scale
    # 将UMAP 的数值 映射到spatial2上
    df_plot['scale_UMAP1'] = scale(df_plot['UMAP1'],
                                          edge_min=df_plot['spatial1'].min(),
                                          edge_max=df_plot['spatial1'].max())
    df_plot['scale_UMAP2'] = scale(df_plot['UMAP2'],
                                          edge_min=df_plot['spatial2'].min(),
                                          edge_max=df_plot['spatial2'].max())
    ## 以0点颠倒
    df_plot['scale_spatial2'] = scale(df_plot['spatial2'].mul(-1),
                                      edge_max=df_plot['spatial2'].max(),
                                      edge_min=df_plot['spatial2'].min())
    # scatter
    ## scatter [spatial]
    [ax.scatter('spatial1', 'scale_spatial2', zs=height, label=label, s=circle_radius,
                marker=marker, c=colormap[label],
                data=df_plot.query("{} == '{}'".format(key, label)),
                **kw_scatter)
        for label in colormap.keys()]
    
    ## scatter [UMAP]
    [ax.scatter('scale_UMAP1', 'scale_UMAP2', zs=0, label=label, s=circle_radius,
                marker=marker, c=colormap[label],
                data=df_plot.query("{} == '{}'".format(key, label)),
                **kw_scatter)
        for label in colormap.keys()]
    
    
    # line
    if not query_3d_line:
        query_3d_line = "{0} == {0}".format(key)
    for i_plot, row_plot in df_plot.query(query_3d_line).iterrows():
        kw_line.update({'color': colormap[row_plot[key]]})
        ax.plot3D(xs=[row_plot['scale_UMAP1'], row_plot['spatial1']],
                  ys=[row_plot['scale_UMAP2'], row_plot['scale_spatial2']],
                  zs=[0, height], alpha=.25, **kw_line)

    # 调节视图
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    xlim, ylim, zlim = ax.get_xlim(), ax.get_ylim(), ax.get_zlim()
    ax.set_axis_off()
    ax.view_init(**kw_view_init)
    
    return df_plot

