#!/usr/bin/env python
# coding: utf-8

# 
# ```bash
# conda activate
# 
# cd ~/link/csMAHN_Spatial
# jupyter nbconvert utils/*.ipynb --to python && rm utils/del_py.py
# jupyter nbconvert utils/scanpy/*.ipynb --to python
# jupyter nbconvert utils/plot/*.ipynb --to python
# jupyter nbconvert init.ipynb --to python
# jupyter nbconvert README.ipynb --to markdown
# 
# 
# :
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


# In[ ]:


# ut.general
Block = ut.Block
json = ut.json

Path = ut.df.Path
np = ut.df.np
pd = ut.df.pd

sc = ut.sc.sc


# In[ ]:


pl = ut.pl.pl
mpl = pl.figure.mpl
plt = pl.figure.plt
sns = pl.figure.sns


# In[4]:


# print(*[_ for _ in  dir(ut)
#         if not _.startswith('__')],
#       sep='\n')

