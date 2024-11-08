#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import torch

sc.settings.verbosity = 0


# In[2]:


import GenKI as gk
from GenKI.preprocesing import build_adata
from GenKI.dataLoader import DataLoader
from GenKI.train import VGAE_trainer
from GenKI import utils

get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[15]:


# subset data as an example

adata = build_adata("../data/adata.h5ad", scale_data=True)
adata
target_cells = 6800
adata.n_obs/target_cells
adata_1 = adata[np.random.choice(adata.obs_names, size=target_cells, replace=True)]


# In[18]:


# load data
target_genes = ['LYAR', 'SSRP1', 'PTPN14', 'TP53', 'YY1']
N = len(target_genes)

for i, gene in enumerate(target_genes[:N]):
    data_wrapper = DataLoader(
                adata_1, # adata object
                target_gene = ["TP53"], # KO gene name
                target_cell = None, # obsname for cell type, if none use all
                obs_label = "ident", # colname for genes
                GRN_file_dir = "GRNs", # folder name for GRNs
                rebuild_GRN = True, # whether build GRN by pcNet
                pcNet_name = "grNet", # GRN file name
                verbose = True, # whether verbose
                n_cpus = 8, # multiprocessing
                )

data_wt = data_wrapper.load_data()
data_ko = data_wrapper.load_kodata()


# In[ ]:


# init trainer

hyperparams = {"epochs": 100, 
               "lr": 7e-4, 
               "beta": 1e-4, 
               "seed": 8096}
log_dir=None 

sensei = VGAE_trainer(data_wt, 
                     epochs=hyperparams["epochs"], 
                     lr=hyperparams["lr"], 
                     log_dir=log_dir, 
                     beta=hyperparams["beta"],
                     seed=hyperparams["seed"],
                     verbose=False,
                     )


# In[ ]:


# %%timeit

sensei.train()


# In[ ]:


# save model

sensei.save_model('adata_genki')


# In[ ]:


# get distance between wt and ko

z_mu_wt, z_std_wt = sensei.get_latent_vars(data_wt)
z_mu_ko, z_std_ko = sensei.get_latent_vars(data_ko)
dis = gk.utils.get_distance(z_mu_ko, z_std_ko, z_mu_wt, z_std_wt, by="KL")
print(dis.shape)


# In[ ]:


# raw ranked gene list

res_raw = utils.get_generank(data_wt, dis, rank=True)
res_raw.head()


# In[ ]:


# if permutation test

null = sensei.pmt(data_ko, n=100, by="KL")
res = utils.get_generank(data_wt, dis, null,)

# save_significant_as = 'gene_list_example')
res

