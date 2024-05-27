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
# echo '[finish] nbconvert'
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


# In[4]:


from utils.general import *


# In[21]:


_temp = np.array([_ for _ in dir(ut.general) if not _.startswith('_')])
print('\n[names import from utils.general]\n'.center(150,'-'))
for _ in np.array_split(_temp,max(2,np.ceil(_temp.size//4))):
    print('  {}'.format(' '.join(_)))
del _temp


# # 包名暴露

# In[ ]:


def module_exists(module_name):
    import importlib.util
    import sys
    return module_name in sys.modules or importlib.util.find_spec(module_name)


# In[ ]:


with Block('[import utils.scanpy]',context={
    'module':'scanpy'.split(',')
}) as context:
    if all([ module_exists(_) for _ in context.module]):
        sc = ut.sc.sc
    else:
        sc = ut.sc


# In[5]:


with Block('[import utils.plot]',context={
    'module':'matplotlib,seaborn'.split(',')
}) as context:
    if all([ module_exists(_) for _ in context.module]):
        pl = ut.pl.pl        
    else:
        pl = ut.pl


# # path

# In[7]:


p_root = Path('~/link/csMAHN_Spatial').expanduser()
p_cache = p_root.joinpath('dataset/cache')

[_.mkdir(parents=True,exist_ok=True) for _ in [p_cache]]

