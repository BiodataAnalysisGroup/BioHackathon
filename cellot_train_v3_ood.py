import sys
from pathlib import Path
from cellot.train.train import train_cellot, train_auto_encoder, train_popalign
from types import SimpleNamespace
from cellot.utils.loaders import load
import scanpy as sc

#from pathlib import Path
#from cellot_train import train_cellot, train_auto_encoder, train_popalign

#mport csv
#from pathlib import Path

#import torch
#import numpy as np
#import random
#import pickle
#from absl import logging
#from absl.flags import FLAGS
#from cellot import losses
#from cellot.utils.loaders import load
#from cellot.models.cellot import compute_loss_f, compute_loss_g, compute_w2_distance
#from cellot.train.summary import Logger
#from cellot.data.utils import cast_loader_to_iterator
#from cellot.models.ae import compute_scgen_shift
#from tqdm import trange





class ConfigNamespace(SimpleNamespace):
    def get(self, key, default=None):
        return getattr(self, key, default)

    def to_dict(self):
        """
        Recursively converts the ConfigNamespace object into a dictionary.
        """
        result = {}
        for key, value in self.__dict__.items():
            if isinstance(value, ConfigNamespace):
                result[key] = value.to_dict()
            else:
                result[key] = value
        return result

    def as_dict(self):
        """Returns the instance as a dictionary for secure use in code"""
        return self.to_dict()

    def __contains__(self, key):
        return key in self.__dict__

# transform a dictionnary in ConfigNamespace
def dict_to_namespace(config_dict):
    return ConfigNamespace(**{k: dict_to_namespace(v) if isinstance(v, dict) else v for k, v in config_dict.items()})

# Function for converting ConfigNamespace objects into a dictionary before use
def convert_to_dict_if_namespace(obj):
    """converting ConfigNamespace objects into a dictionary"""
    if isinstance(obj, ConfigNamespace):
        return obj.to_dict()
    elif isinstance(obj, dict):
        return {k: convert_to_dict_if_namespace(v) for k, v in obj.items()}
    else:
        return obj



def run_cellot_training(task_config, model_config, train_type="cellot", outdir='./output'):
    config = {
        'training': {
            'n_iters': task_config.get('epochs', 10),
            'logs_freq': 10,
            'eval_freq': 20,
            'cache_freq': 100,
            'n_inner_iters': 1
        },
        'data': {
            'condition': task_config.get('condition', 'drug'),
            'source': task_config.get('source', 'control'),
            'target': task_config.get('target', 'stim'),
            'type': task_config.get('type', 'cell'),
            'path': task_config.get('dataset', '')
        },
        'model': model_config,
        'datasplit': {
            'groupby': task_config.get('datasplit_groupby', 'condition'),
            'name': task_config.get('datasplit_name', 'train_test'),
            'test_size': task_config.get('datasplit_test_size', 0.2),
            'random_state': task_config.get('datasplit_random_state', 0),
            'holdout': task_config.get('datasplit_holdout', None),
            'key': task_config.get('key', None),
            'mode': task_config.get('datasplit_mode', 'iid'),
            'subset': None
        },
        'dataloader': {
            'batch_size': task_config.get('batch_size', 64),
            'shuffle': task_config.get('shuffle', True)
        }
    }

    # Transform in a ConfigNamespace
    config_ns = dict_to_namespace(config)

    # Convert outdir to Path object to ensure compatibility
    outdir = Path(outdir)

    # Call the model function
    if train_type == 'cellot':
        train_cellot(outdir, config_ns)
    elif train_type == 'auto_encoder':
        train_auto_encoder(outdir, config_ns)
    elif train_type == 'popalign':
        train_popalign(outdir, config_ns)
    else:
        raise ValueError("Train type not supported: {}".format(train_type))

    print(f"Training complete for {train_type} model at {outdir}.")


import os
import anndata
#from cellot_train_v2 import run_cellot_training  # Ensure this is the path to the training function

# Load the AnnData dataset
adata = anndata.read_h5ad("C:\\Users\\Shadow\\Desktop\\BioHack24\\scPRAM\\processed_datasets_all\\datasets\\scrna-lupuspatients\\kang-hvg.h5ad")  # Replace with the actual path

# Extract unique cell types
cell_types = adata.obs['cell_type'].unique()

# Directory to save the models
output_dir = '.\\output_ood_models'


for cell_type in cell_types:
    # Set up holdout as the current cell type for OOD
    holdout = cell_type
    
    # Define the output directory for the model
    model_dir = os.path.join(output_dir, f"{cell_type}_ood")
    os.makedirs(model_dir, exist_ok=True)
    
    # Define the task configuration for training
    task_config = {
        'dataset': "C:\\Users\\Shadow\\Desktop\\BioHack24\\scPRAM\\processed_datasets_all\\datasets\\scrna-lupuspatients\\kang-hvg.h5ad",
        'condition': 'condition',      # Column defining the conditions
        'source': 'ctrl',              # Control condition as source
        'target': 'stim',              # Stimulated condition as target
        'type': 'cell',                # Specifies type, assuming it's for cells
        'epochs': 100000,              # Set the number of epochs
        'batch_size': 128,             # Batch size
        'shuffle': True,               # Shuffle the data
        'datasplit_groupby': ['cell_type','condition'],  # Grouping by 'condition' for OOD
        'datasplit_name': 'toggle_ood',
        'key' : 'cell_type',
        'datasplit_mode': 'ood',           # Set mode to 'ood'
        'datasplit_holdout': holdout,      # Specify holdout cell type
        'datasplit_test_size': 0.3,        # Test split size
        'datasplit_random_state': 0        # Random seed for reproducibility
    }

    # Define the model configuration
    model_config = {
        'input_dim': 1000,
        'name': 'cellot',
        'hidden_units': [64, 64, 64, 64],
        'latent_dim': 100,
        'softplus_W_kernels': False,
        'g': {  
            'fnorm_penalty': 1
        },
        'kernel_init_fxn': {
            'b': 0.1,
            'name': 'uniform'
        },
        'optim': {
            'optimizer': 'Adam',
            'lr': 0.0001,
            'beta1': 0.5,
            'beta2': 0.9,
            'weight_decay': 0
        },
        'training': {
            'n_iters': 100000,       # Total number of iterations
            'n_inner_iters': 1,      # Number of inner iterations
            'cache_freq': 50,       # Frequency to cache model state
            'eval_freq': 20,         # Frequency for evaluations
            'logs_freq': 10          # Logging frequency
        }
    }
    
    print(f"Training model for holdout cell_type: {cell_type}")
    
    # Run the training with specified configurations and train_type='cellot'
    run_cellot_training(
        task_config=task_config,
        model_config=model_config,
        train_type='cellot',
        outdir=model_dir  # Save output to the model-specific directory
    )
    
    print(f"Completed training for holdout cell_type: {cell_type}")



#exemple


test_data_path = "C:\\Users\\Shadow\\Desktop\\BioHack24\\scPRAM\\processed_datasets_all\\datasets\\scrna-lupuspatients\\kang-hvg.h5ad"
adata = sc.read_h5ad(test_data_path)
input_dim = adata.shape[1]
print("Number of variables (input_dim) :", input_dim)





import pandas as pd
import matplotlib.pyplot as plt

# load losses data
loss_data = pd.read_csv('.\\output_ood_models\\loss_tracking.csv')

# Check for column
if all(col in loss_data.columns for col in ['Step', 'Loss_G', 'Loss_F']):
    # plot
    plt.figure(figsize=(10, 6))
    plt.plot(loss_data['Step'], loss_data['Loss_G'], label='Loss_G')
    plt.plot(loss_data['Step'], loss_data['Loss_F'], label='Loss_F')
    plt.xlabel('Step')
    plt.ylabel('Loss')
    plt.yscale('log')
    plt.legend()
    plt.title('Evolution of Losses (Loss G and Loss F) During Training')
    plt.show()
else:
    print("No 'Step', 'Loss_G', ou 'Loss_F' in the CSV file.")

