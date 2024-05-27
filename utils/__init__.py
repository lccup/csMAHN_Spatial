#!/usr/bin/env python
# coding: utf-8

# ```bash
# 
# conda activate
# 
# cd ~/link/csMAHN_Spatial/utils
# {
# jupyter nbconvert *.ipynb --to python && rm del_py.py
# jupyter nbconvert scanpy/*.ipynb --to python
# jupyter nbconvert plot/*.ipynb --to python
# clear
# echo '----------------------------------------'
# echo '[finish] nbconvert'
# echo '----------------------------------------'
# }
# ```
# 
# 
# ```mermaid
# graph LR
#     general[general]
#     general --> df[df]
#     general --> arr[arr]
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


# In[2]:


def module_exists(module_name):
    import importlib.util
    import sys
    return module_name in sys.modules or importlib.util.find_spec(module_name)


# In[3]:


from utils import general
from utils.general import *


# In[ ]:


from utils import arr
from utils import df


# In[ ]:


with Block('[import utils.scanpy]',context={
    'module':'scanpy'.split(',')
}) as context:
    if all([ module_exists(_) for _ in context.module]):
        from utils import scanpy as sc
    else:
        sc = '[module has not installed] {}'.format(','.join(context.module))


# In[ ]:


with Block('[import utils.plot]',context={
    'module':'matplotlib,seaborn'.split(',')
}) as context:
    if all([ module_exists(_) for _ in context.module]):
        import utils.plot as pl
    else:
        sc = '[module has not installed] {}'.format(','.join(context.module))

