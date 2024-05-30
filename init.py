#!/usr/bin/env python
# coding: utf-8

# 
# ```bash
# conda activate
# 
# cd ~/link/csMAHN_Spatial
# {
# jupyter nbconvert utils/*.ipynb --to python && rm utils/del_py.py
# jupyter nbconvert utils/__init__.ipynb --to html
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


from IPython.display import display

import utils as ut
from utils.general import *
# print(ut.__doc__)


# In[2]:


_temp = np.array([_ for _ in dir(ut.general) if not _.startswith('_')])
print('\n[names import from utils.general]\n'.center(150,'-'))
for _ in np.array_split(_temp,max(2,np.ceil(_temp.size//4))):
    print('  {}'.format(' '.join(_)))
del _temp


# # 包名暴露

# In[3]:


with Block('[import utils.scanpy]',context={
    'module':'scanpy'.split(',')
}) as context:
    if all([ module_exists(_) for _ in context.module]):
        sc = ut.sc.sc
    else:
        sc = ut.sc


# In[4]:


with Block('[import utils.plot]',context={
    'module':'matplotlib,seaborn'.split(',')
}) as context:
    if all([ module_exists(_) for _ in context.module]):
        pl = ut.pl   
    else:
        pl = ut.pl


# # path

# In[11]:


p_root = Path('~/link/csMAHN_Spatial').expanduser()
p_cache = p_root.joinpath('dataset/cache')
p_res = p_root.joinpath('res')

[_.mkdir(parents=True,exist_ok=True) for _ in [p_cache,p_res]]


# In[6]:


_df_species_name = pd.DataFrame({
'simple':'hs,ma,mm,dr'.split(','),
'common':'human,macaque,mouse,zebrafish'.split(','),
'official':'Homo sapiens,Macaca,Mus musculus,Danio rerio'.split(',')
})
def convert_species_name(value_query,
                         key_query='simple',
                         key_return='common',
                         df_species_name=_df_species_name):
    assert key_query in df_species_name.columns
    assert key_return in df_species_name.columns
    temp = df_species_name.query("`{}` == '{}'".format(key_query,value_query))
    assert temp.shape[0] == 1,'[Error] can not get unique item'
    return temp[key_return].to_numpy()[0]


# # res_publish func

# In[7]:


import sys
from pathlib import Path
p_root = Path('~/link/res_publish').expanduser()
sys.path.append(str(p_root)) if not str(p_root) in sys.path else None


# In[8]:


import func


# In[9]:


func.run_cross_species_models

