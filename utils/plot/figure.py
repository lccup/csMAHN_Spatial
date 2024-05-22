#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""figure

绘图初始化
分图

> variable
    rc_frame
> function
    subplots_get_fig_axs
"""


# In[2]:


from  pathlib import Path

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


# In[3]:


rc_default = {}
# fig
rc_default.update({
    'figure.dpi':300
})
# font
rc_default.update({
    'font.family':'arial',
    'font.size':6  # 磅（points）
})
# axes
rc_default.update({
    'axes.facecolor':'white',
    'axes.labelsize':rc_default['font.size'],
    'axes.titlesize':rc_default['font.size'] + 2,
    'axes.titleweight':'bold',
    # axes.edge
    'axes.edgecolor': 'white',
    'axes.linewidth':.5,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'axes.spines.left': False,
    'axes.spines.bottom': False,
})
# legend
rc_default.update({
    'legend.title_fontsize':rc_default['font.size']
})
# tick
rc_default.update({
    'xtick.major.size': 0,
    'xtick.major.width':.2,
    'xtick.color': 'black',
    'xtick.major.pad':1,
    'xtick.labelsize':rc_default['font.size'],

    'ytick.major.size': 0,
    'ytick.major.width':.2,
    'ytick.color': 'black',
    'ytick.major.pad':1,
    'ytick.labelsize':rc_default['font.size']
})
# patch
rc_default.update({
    'patch.edgecolor':'white',
    'patch.linewidth':.5
})


# In[9]:


from  pathlib import Path
print(Path('.').absolute())


# In[4]:


mpl.font_manager.fontManager.addfont(
    Path(__file__).parent.joinpath('font/arial.ttf'))
sns.set_theme(style="ticks",rc=mpl.rcParams)
mpl.rcParams.update(rc_default)
rc_frame = rc_default.copy()
rc_frame.update({
    'axes.edgecolor': 'black',
    'axes.linewidth':.5,
    "axes.spines.left": True,
    "axes.spines.bottom": True,
    'xtick.major.size':2,
    'ytick.major.size':2
})


# In[5]:


def subplots_get_fig_axs(nrows=1, ncols=1,
                             ratio_nrows=2, ratio_ncols=2):
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(
        ratio_nrows*ncols, ratio_ncols*nrows)
    )
    axs = np.ravel(axs)
    if axs.size == 1:
        axs = axs[0]

    return fig, axs

