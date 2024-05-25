#!/usr/bin/env python
# coding: utf-8

# 
# ```bash
# conda activate
# 
# cd ~/link/csMAHN_Spatial
# {
# jupyter nbconvert utils/*.ipynb --to python && rm utils/del_py.py
# jupyter nbconvert utils/scanpy/*.ipynb --to python
# jupyter nbconvert utils/plot/*.ipynb --to python
# jupyter nbconvert init.ipynb --to python
# jupyter nbconvert README.ipynb --to markdown
# clear
# echo '----------------------------------------'
# echo '[finish]'
# echo '----------------------------------------'
# }
# 
# ```

# In[1]:


from pathlib import Path
Path('/public/workspace/licanchengup/link/').parent


# In[2]:


# utils中有，就不再导入了，太耗时间了
# from pathlib import Path
# import numpy as np
# import pandas as pd

import utils as ut
print(ut.__doc__)


# In[3]:


from utils.general import *


# # 包名暴露

# In[4]:


Path = ut.df.Path
np = ut.df.np
pd = ut.df.pd

mpl = ut.pl.figure.mpl
plt = ut.pl.figure.plt
sns = ut.pl.figure.sns

sc = ut.sc.sc


# In[5]:


pl = ut.pl.pl


# In[6]:


# print(*[_ for _ in  dir(ut)
#         if not _.startswith('__')],
#       sep='\n')


# # path

# In[7]:


p_root = Path('~/link/csMAHN_Spatial').expanduser()
p_cache = p_root.joinpath('dataset/cache')

[_.mkdir(parents=True,exist_ok=True) for _ in [p_cache]]

