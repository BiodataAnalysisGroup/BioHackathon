import torch
from pathlib import Path
import pandas as pd
import scanpy as sc
from cellot.models.cellot import load_cellot_model, compute_loss_g, compute_loss_f, compute_w2_distance
import anndata
import numpy as np
import os
import matplotlib.pyplot as plt
from types import SimpleNamespace
from sklearn.metrics import r2_score
print(os.getcwd())



class ConfigNamespace(SimpleNamespace):
    def get(self, key, default=None):
        return getattr(self, key, default)

    def to_dict(self):
        """
        Convertit récursivement l'objet ConfigNamespace en dictionnaire.
        """
        result = {}
        for key, value in self.__dict__.items():
            if isinstance(value, ConfigNamespace):
                result[key] = value.to_dict()
            else:
                result[key] = value
        return result

    def as_dict(self):
        """Renvoie l'instance sous forme de dictionnaire pour une utilisation sécurisée dans le code."""
        return self.to_dict()

    def __contains__(self, key):
        return key in self.__dict__

# Fonction utilitaire pour transformer un dictionnaire en ConfigNamespace
def dict_to_namespace(config_dict):
    return ConfigNamespace(**{k: dict_to_namespace(v) if isinstance(v, dict) else v for k, v in config_dict.items()})

# Fonction pour convertir les objets ConfigNamespace en dictionnaire avant l'utilisation
def convert_to_dict_if_namespace(obj):
    """Convertit un ConfigNamespace en dictionnaire, récursivement si nécessaire."""
    if isinstance(obj, ConfigNamespace):
        return obj.to_dict()
    elif isinstance(obj, dict):
        return {k: convert_to_dict_if_namespace(v) for k, v in obj.items()}
    else:
        return obj



def load_test_data(test_data_path, config):
    test_data = sc.read(test_data_path)
    source_data = test_data[test_data.obs[config.data.condition] == config.data.source]
    target_data = test_data[test_data.obs[config.data.condition] == config.data.target]
    source_tensor = torch.tensor(source_data.X.toarray(), dtype=torch.float32, requires_grad=True)
    target_tensor = torch.tensor(target_data.X.toarray(), dtype=torch.float32)
    return list(zip(source_tensor, target_tensor))


def create_anndata_with_predictions(config, model_path, original_data):
    # Load the model
    checkpoint = torch.load(model_path)
    (f, g), _ = load_cellot_model(config)
    f.load_state_dict(checkpoint['f_state'])
    g.load_state_dict(checkpoint['g_state'])

    # Set the model to evaluation mode
    f.eval()
    g.eval()

    # Filter for source cells (ctrl condition)
    source_data = original_data[original_data.obs[config.data.condition] == config.data.source]
    
    # Convert source data to tensor and set requires_grad for all tensors
    source_tensor = torch.tensor(
        source_data.X.toarray() if hasattr(source_data.X, "toarray") else source_data.X, 
        dtype=torch.float32, 
        requires_grad=True
    )

    # Step 1: Verify for NaNs in the source_tensor
    print(f"Step 1: Nombre de NaN dans source_tensor : {torch.isnan(source_tensor).sum().item()}")

    # Store predicted cells
    predicted_cells = []

    # Step 2: Obtain predictions for each source cell and check for NaNs in the prediction
    for i, source in enumerate(source_tensor):
        with torch.set_grad_enabled(True):  # Ensure grad tracking is enabled
            source = source.unsqueeze(0)  # Ensure correct shape
            predicted = g.transport(source)  # Transport function
            
            # Check if prediction contains NaNs
            if torch.isnan(predicted).any():
                print(f"Step 2: NaNs detected in prediction for cell {i}")
            else:
                print(f"Step 2: Prediction successful for cell {i}")

            predicted_cells.append(predicted.detach().numpy())  # Detach after prediction

    # Stack predictions into an array
    predicted_data_matrix = np.vstack(predicted_cells)
    # Normalize predicted data to match original data's scale
    # This assumes that the original data has been normalized with scanpy's `normalize_total`
    predicted_adata = anndata.AnnData(X=predicted_data_matrix)
    #sc.pp.normalize_total(predicted_adata, target_sum=1e4)
    #sc.pp.log1p(predicted_adata)
    #predicted_data_matrix = predicted_adata.X


    # Step 3: Verify if predicted_data_matrix contains NaNs after prediction loop
    print(f"Step 3: Nombre de NaN dans predicted_data_matrix : {np.isnan(predicted_data_matrix).sum()}")

    # Optional: Replace NaNs in `predicted_data_matrix` if needed
    # predicted_data_matrix = np.nan_to_num(predicted_data_matrix)
    # print(f"Step 3b: NaNs replaced. Nombre de NaN dans predicted_data_matrix après remplacement : {np.isnan(predicted_data_matrix).sum()}")

    # Combine predicted data matrix with the existing data matrix, ensuring no duplication with 'ctrl'
    original_data_matrix = (
        original_data.X.toarray() if hasattr(original_data.X, "toarray") else original_data.X
    )
    
    # Step 4: Check for NaNs in the original data matrix
    print(f"Step 4: Nombre de NaN dans original_data_matrix : {np.isnan(original_data_matrix).sum()}")

    # Step 5: Combine matrices and check for NaNs in the combined data
    combined_data = np.vstack([original_data_matrix, predicted_data_matrix])
    print(f"Step 5: Nombre de NaN dans combined_data après concaténation : {np.isnan(combined_data).sum()}")

    # Copy original metadata and create labels for predictions
    combined_obs = original_data.obs.copy()

    # Generate a new observation dataframe for predicted cells based on source cells but labeled as 'predicted'
    predicted_obs = source_data.obs.copy()
    predicted_obs[config.data.condition] = 'predicted'  # Set new condition
    predicted_obs.index = [f"pred_{i}" for i in range(len(predicted_cells))]  # Unique indices for predictions

    # Concatenate the original observations with the newly created predicted observations
    combined_obs = pd.concat([combined_obs, predicted_obs])

    # Final AnnData object with original and predicted cells
    anndata_with_predictions = anndata.AnnData(
        X=combined_data,
        obs=combined_obs,
        var=original_data.var
    )


    # Ensure observation names are unique
    anndata_with_predictions.obs_names_make_unique()

    # Step 6: Check if the final AnnData object contains NaNs in X
    print(f"Step 6: Nombre de NaN dans anndata_with_predictions.X : {np.isnan(anndata_with_predictions.X).sum()}")

    # Optional: If desired, set raw attribute for the AnnData object
    anndata_with_predictions.raw = anndata_with_predictions.copy()

    return anndata_with_predictions


# Load the dataset
dataset_path = "C:\\Users\\Shadow\\Desktop\\BioHack24\\scPRAM\\processed_datasets_all\\datasets\\scrna-lupuspatients\\kang-hvg.h5ad"
adata = sc.read_h5ad(dataset_path)
cell_types = adata.obs['cell_type'].unique()
output_dir = ".\\output_ood_models"

for cell_type in cell_types:
    model_dir = os.path.join(output_dir, f"{cell_type}_ood")
    model_path = Path(model_dir) / "cache" / "model.pt"
    holdout = cell_type
    
    # Define task configuration for evaluation
    task_config = {
        'dataset': dataset_path,
        'condition': 'condition',
        'source': 'ctrl',
        'target': 'stim',
        'type': 'cell',
        'batch_size': 128,
        'shuffle': True,
        'datasplit_groupby': ['cell_type','condition'],
        'datasplit_name': 'toggle_ood',
        'key' : 'cell_type',# Name for this data split
        'datasplit_mode': 'ood',           # Set mode to 'ood'
        'datasplit_holdout': holdout,      # Specify holdout cell type
        'datasplit_test_size': 0.3,
        'datasplit_random_state': 0
    }

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
            'n_iters': 100000,
            'n_inner_iters': 1,
            'cache_freq': 50,
            'eval_freq': 20,
            'logs_freq': 10
        }
    }

    config = {
        'training': model_config['training'],
        'data': task_config,
        'model': model_config,
        'datasplit': {
            'groupby': task_config['datasplit_groupby'],
            'name': task_config['datasplit_name'],
            'test_size': task_config['datasplit_test_size'],
            'random_state': task_config['datasplit_random_state'],
            'holdout': task_config.get('datasplit_holdout', None),
            'key': task_config.get('key', None),
            'mode': task_config.get('datasplit_mode', 'iid'),            
            'subset': None
        },
        'dataloader': {
            'batch_size': task_config['batch_size'],
            'shuffle': task_config['shuffle']
        }
    }

    config_ns = dict_to_namespace(config)
    
    # Evaluate model and create AnnData with predictions
    anndata_with_predictions = create_anndata_with_predictions(config_ns, model_path, adata)
 
    
    # Filter only for the cells of the holdout type and their predictions
    holdout_cells = anndata_with_predictions[
        (anndata_with_predictions.obs['cell_type'] == cell_type) 
    ]

    # Visualization with PCA
    print(f"Evaluating PCA and UMAP for holdout cell_type: {cell_type}")
    sc.tl.pca(holdout_cells, svd_solver="arpack")
    sc.pl.pca(holdout_cells, color="condition", title=f"PCA for {cell_type}")

    # Visualization with UMAP
    sc.pp.neighbors(holdout_cells)
    sc.tl.umap(holdout_cells)
    sc.pl.umap(holdout_cells, color="condition", title=f"UMAP for {cell_type}") 
    



