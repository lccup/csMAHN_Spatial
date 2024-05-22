#!/usr/bin/env python
# coding: utf-8

# In[48]:


"""colormap
颜色

> function
    get
    show
"""


# In[2]:


from pathlib import Path
import json

from utils.plot.figure import *


# In[3]:


palettes = json.loads(Path(__file__).parent\
                      .joinpath('scanpy.plotting.palettes.json').read_text())
palettes.update({k: v.split(',') for k, v in palettes.items()
                 if k.startswith('default')})


# In[ ]:


def get_color(serise):
    palette = None
    serise = pd.Series(pd.Series(serise).unique())
    if serise.size <= 20:
        palette = palettes['default_20']
    elif serise.size <= 28:
        palette = palettes['default_28']
    elif serise.size <= len(palettes['default_102']):  # 103 colors
        palette = palettes['default_102']
    else:
        raise Exception("[categories too long] {}".format(serise.size))
    return palette[:serise.size]

def get(serise, color_missing_value="lightgray",
        offset=2, filter_offset=True):

    serise = pd.Series(pd.Series(serise).unique())
    has_missing_value = serise.isna().any()
    serise = pd.Series(np.concatenate(
        (['_{}'.format(i) for i in range(offset)], serise.dropna().astype(str))))

    palette = get_color(serise)

    colormap = {k: v for k, v in zip(serise, palette)}
    if has_missing_value:
        colormap.update({'nan': color_missing_value})

    if filter_offset:
        colormap = {
            k: v
            for _, (k, v) in zip(
                ~pd.Series(colormap.keys()).str.match('_\\d+'),
                colormap.items())
            if _
        }
    return colormap

def show(color_map,
         marker='.', size=40, fontdict={
             'horizontalalignment': 'left',
             'verticalalignment': 'center'},
         ax=None, return_fig=False):

    if ax:
        fig = ax.figure
    else:
        fig, ax = subplots_get_fig_axs()
    if isinstance(marker, str):
        marker = np.repeat(marker, len(color_map.keys()))
    for i, ((k, v), m) in enumerate(zip(color_map.items(), marker)):
        ax.scatter(0, len(color_map.keys())-i,
                   label=k, c=v, s=size, marker=m)
        ax.text(0.1, len(color_map.keys())-i, k, fontdict=fontdict)
    ax.set_xlim(-0.05, 0.25)
    ax.set_ymargin(.5)
    ax.set_axis_off()
    if return_fig:
        return fig


# In[ ]:


if __name__ == '__main__':
    fig,axs = subplots_get_fig_axs(ncols=3)
    for i in range(3):
        show(get([' '],offset=i+3,filter_offset=False),ax=axs[i])

