#!/usr/bin/env python
# coding: utf-8

# ```bash
# 
# conda activate
# 
# cd ~/link/csMAHN_Spatial
# jupyter nbconvert utils/*.ipynb --to python && rm utils/del_py.py
# jupyter nbconvert utils/scanpy/*.ipynb --to python
# jupyter nbconvert utils/plot/*.ipynb --to python
# :
# ```
# 
# 
# ```mermaid
# graph LR
#     general[general]
#     df[df]
#     general --> df
# ```

# In[ ]:


"""
LCC 的python 函数库
为实验室留下些什么吧
0.0.1 2024年5月22日10:48:05

import utils as ut
# help(ut)
print(ut.__doc__)
print(ut.__version__)
"""


# In[ ]:


__version__ = '0.0.1'


# In[ ]:


from utils.general import *


# In[ ]:


from utils import arr
from utils import df


# In[ ]:


from utils import scanpy as sc


# In[ ]:


import utils.plot as pl

