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
> class
    figure_A4Page
"""


# In[2]:


from utils.general import Path,np,pd
from utils.general import mpl,plt,sns


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

# legend
rc_default.update({
    'legend.fontsize':rc_default['font.size'],
    'legend.frameon':False,
    'legend.markerscale':20
})


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
                         ratio_nrows=2, ratio_ncols=2,
                         width_ratios=None, height_ratios=None,
                         kw_subplot=None, kw_gridspec=None,
                         kw_fig={}, ravel=True, kw_ravel={}):
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols,
                            width_ratios=width_ratios, height_ratios=height_ratios,
                            gridspec_kw=kw_gridspec, subplot_kw=kw_subplot,
                            figsize=(ratio_nrows*ncols, ratio_ncols*nrows), **kw_fig)

    if ravel:
        axs = np.ravel(axs, **kw_ravel)
    if axs.size == 1:
        axs = axs[0]

    return fig, axs


# # figure_A4Page
# 
# ```python
# # A4大小的figure,使用网格控制子图
# display([8.27*.97, 11.96*.96])
# display([8.27*.97/.25, 11.96*.96/.25])
# # 以0.25为单位1，将width分为32等分，higth分为46等分
# fig = plt.figure(figsize=(8.27, 11.69))
# spec = fig.add_gridspec(nrows=46, ncols=32,
#                     left=0.03, right=1,  # 设置边距
#                     bottom=0.02, top=0.98,  # 设置边距
#                     wspace=0, hspace=0)  # 设置子图间距
# ```

# In[ ]:


class figure_A4Page:
    """
A4大小的figure,使用网格控制子图
display([8.27*.97, 11.96*.96])
display([8.27*.97/.25, 11.96*.96/.25])
# 以0.25为单位1，将width分为32等分，higth分为46等分
fig = plt.figure(figsize=(8.27, 11.69))
spec = fig.add_gridspec(nrows=46, ncols=32,
                    left=0.03, right=1,  # 设置边距
                    bottom=0.02, top=0.98,  # 设置边距
                    wspace=0, hspace=0)  # 设置子图间距
    """
    def __init__(self):
        fig = plt.figure(figsize=(8.27, 11.69))
        spec = fig.add_gridspec(nrows=46, ncols=32,
                            left=0.03, right=1,  # 设置边距
                            bottom=0.02, top=0.98,  # 设置边距
                            wspace=0, hspace=0)  # 设置子图间距
        self.fig = fig
        self.spec = spec

    def add_grid(self,alpha=0.15):
        """添加网格
        """

        for i,(v) in enumerate(np.linspace(0.03,1,32+1)):
            # |
            self.fig.add_artist(mpl.lines.Line2D([v,v],[0.02,0.98],linestyle='--',alpha=alpha,c = '#00BFFF'))
            if i >= 32 or i%2 ==1:
                continue
            for y in [0.1,0.5,0.9]:
                self.fig.text(v+0.01,y,'x={:.0f}'.format(i) if i==0 or i== 16 else '{:.0f}'.format(i),
                         fontdict={'fontsize':10,'color':'blue','alpha':alpha})
        fontdict={'fontsize':10,'color':'blue','alpha':alpha}
        for i,(v) in enumerate(np.linspace(0.98,0.02,46 + 1)):
            # ------
            self.fig.add_artist(mpl.lines.Line2D([0.03,1],[v,v],linestyle='--',alpha=alpha,c = '#D02090'))
            if i >= 46 or i%2 ==1:
                continue
            for x in [0.1,0.5,0.9]:
                self.fig.text(x,v-0.01,'y={:.0f}'.format(i) if i==0 or i== 24 else '{:.0f}'.format(i),
                         fontdict={'fontsize':10,'color':'red','alpha':alpha})

    def add_grid_absolute_coordinate(self, alpha=0.15):
        """标出self.fig的绝对坐标系
    """
        for i, (v) in enumerate(np.linspace(0, 1, 10+1)):
            # |
            self.fig.add_artist(
                mpl.lines.Line2D([v, v], [0, 1], linestyle='--', alpha=alpha, c='#2F4F4F'))
            for y in [0.1, 0.5, 0.9]:
                self.fig.text(
                    v+0.01,
                    y,
                    'x=.{:.0f}'.format(i) if v == 0 or v == .5 else '.{:.0f}'.format(i),
                    fontdict={
                        'fontsize': 10,
                        'color': 'black',
                        'alpha': alpha})
        fontdict = {'fontsize': 10, 'color': 'blue', 'alpha': alpha}
    
        for i, (v) in enumerate(np.linspace(0, 1, 10 + 1)):
            # ------
            self.fig.add_artist(
                mpl.lines.Line2D([0, 1], [v, v], linestyle='--', alpha=alpha, c='#FFA500'))
            for x in [0.1, 0.5, 0.9]:
                self.fig.text(
                    x,
                    v+0.01,
                    'y=.{:.0f}'.format(i) if v == 0 or v == .5 else '.{:.0f}'.format(i),
                    fontdict={
                        'fontsize': 10,
                        'color': 'orange',
                        'alpha': alpha})

    def subset_spec(self, x, y, x_offest=1, y_offest=1):
        return self.spec[y:y+y_offest, x:x+x_offest]

    def add_ax_with_spec(self, x, y, x_offest=1, y_offest=1,
                        *args, **kwargs):
        return self.fig.add_subplot(self.subset_spec(x, y,
            x_offest, y_offest),*args, **kwargs)
    
    def add_text_with_ax(self,ax,text,x=0,y=1,fontdict={
            'fontsize': 14,'fontweight': 'bold'}):
        ax.text(x, y, text, fontdict=fontdict)
        ax.set_axis_off()
    def save_as_pdf(self,p,close=True,**kwargs):
        from matplotlib.backends.backend_pdf import PdfPages
        with PdfPages(p) as pdf:
            pdf.savefig(self.fig,**kwargs)
            plt.close('all') if close else None

