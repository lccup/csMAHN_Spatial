#!/usr/bin/env python
# coding: utf-8

# ```bash
# conda activate
# 
# cd ~/link/csMAHN_Spatial
# jupyter nbconvert utils/plot/*.ipynb --to python
# 
# :
# ```
# 
# ```mermaid
# graph LR;
#     __init__{{__init__}};
#     general[general];
#     general[general] -.-> __init__;
# 
# 
#     subgraph plot
#         plot.__init__{{__init__}};
#         plot.figure[figure];
#         plot.colormap[colormap];
#         plot.pl[pl];
#         plot.figure --> plot.colormap
#         plot.figure -.-> plot.pl
#         plot.colormap -.-> plot.pl
#         plot.pl -.-> plot.__init__
#     end
#     general --> plot.figure
#     plot.__init__ -.-> __init__
# 
# ```

# In[2]:


"""utils.plot
借由matplotlib和seaborn实现
"""


# In[ ]:


from utils.plot.pl import *

