#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from pathlib import Path
import json as json

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


# In[ ]:


rng = np.random.default_rng()


# In[ ]:


def show_dict_key(data, tag='', sort_key=True):
    print("\n>{}['']".format(tag).ljust(75, '-'))
    ks = list(data.keys())
    if sort_key:
        ks = np.sort(ks)
    print(*['\t{}'.format(k) for k in ks], sep='\n')


# # Block
# 
# `块` 用于将代码分块
# 
# 因为在notebook中用注释将代码分块，无法实现代码的折叠，使得代码较为混乱
# 
# 故通过`with Block():`组合构造出可以折叠的with语句块，进而提高代码的可读性
# 
# + `with Block():`内部并未与外部进行隔离
# 
# + 实现了计时功能
# + 实现了上下文功能

# In[ ]:


class Block:
    """用于在notebook中给代码划分区域(块),从而使代码能够折叠

# 上下文功能
with Block('context',context={
    'a':'Block','b':':','c':'Hello World'
}) as context:
    print('inner')
    print('\t',' '.join(context.context.values()))
    print('\t',context.a,context.b,context.c)
# output
## inner
## 	 Block : Hello World
## 	 Block : Hello World

# 计时功能
import time
with Block('test',show_comment=True):
    print('inner')
    time.sleep(2)
# output
## [start][test] 0509-00:20:47------------------------------------------------
## inner
## [end][test] 0509-00:20:49--------------------------------------------------
        2.002 s used
    """

    def __init__(self, comment, show_comment=False,context={}):
        self.comment = comment
        self.show_comment = show_comment
        self.context = context
        self.time = 0

    def __enter__(self):
        if self.show_comment:
            self.time = time.time()
            print("[start][{}] {}".format(
                self.comment,
                time.strftime('%m%d-%H:%M:%S',
                              time.localtime())).ljust(75, '-'))
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.show_comment:
            print(
                "[end][{}] {}".format(
                    self.comment,
                    time.strftime(
                        '%m%d-%H:%M:%S',
                        time.localtime())).ljust(
                    75,
                    '-'))
            print("\t{:.3f} s used".format(time.time()-self.time))
        # 释放content
        self.context = None
    def __str__(self):
        return """Block
\tcomment     : {}
\tcontext_key : {}
""".format(self.comment,','.join(self.context.keys()))
    
    # 对类及其实例未定义的属性有效
    # 若name 不存在于 self.__dir__中,则调用__getattr__
    def __getattr__(self,name):
        cls = type(self)
        res = self.context.setdefault(name,None)
        if res:
            return res
        else:
            raise AttributeError(
                '{.__name__!r} object has no attribute {!r}'\
                .format(cls, name))

