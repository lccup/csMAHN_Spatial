#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""df
用于操作pandas.DataFrame
> function
    gourp_add
    iterdir
    apply_merge_field
"""


# In[ ]:


from utils.general import Path,np,pd


# In[ ]:


def show(df,n=2):
    from IPython.display import display
    display(df.head(n),df.shape)


def iterdir(p,path_match='',path_match_filter=[],select='f'):
    p = Path(p)
    assert p.exists(),'[not exists] {}'.format(p)
    assert p.is_dir(),'[Error] p is not a dir'
    
    res = pd.DataFrame({'path':p.iterdir()})
    
    # select
    if select == 'file' or select[0] ==  'f':
        res = res[res['path'].apply(lambda x:x.is_file())]
    elif select == 'dir' or select[0] ==  'd':
        res = res[res['path'].apply(lambda x:x.is_dir())]
    else:
        # file and dir
        pass
    assert res.shape[0] > 0,'[Error] no item'
    
    if path_match:
        res = res[res['path'].apply(lambda x:x.match(path_match))]
 
    if path_match_filter and isinstance(path_match_filter,str):
        path_match_filter = [path_match_filter]
    for _ in path_match_filter:
        res = res[res['path'].apply(lambda x:not x.match(_))]
    
    res['name'] = res['path'].apply(lambda x:x.name)
    res.index = np.arange(res.shape[0])
    return res.copy()

def group_agg(df,groupby_list,agg_dict=None,dropna=True,
        reindex=True,rename_dict=None):
    if None is agg_dict:
        agg_dict = {groupby_list[-1]: ['count']}
    res = df.groupby(
        groupby_list,
        dropna=dropna,
        observed=False).agg(agg_dict)
    if reindex:
        res.columns = ["_".join(i) for i in res.columns]
        res = res.index.to_frame().join(res)
        res.index = np.arange(res.shape[0])
    if isinstance(rename_dict, dict):
        res = res.rename(columns=lambda k: rename_dict.setdefault(k, k))
    return res

def apply_merge_field(df, str_format):
    return df.apply(lambda row: str_format.format(**row), axis=1)


