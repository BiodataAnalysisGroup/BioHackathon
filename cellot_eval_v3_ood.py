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
    # Normalize predicted data to match original data's scale ?
    # This assumes that the original data has been normalized with scanpy's `normalize_total`
    predicted_adata = anndata.AnnData(X=predicted_data_matrix)



    # Step 3: Verify if predicted_data_matrix contains NaNs after prediction loop
    print(f"Step 3: Nombre de NaN dans predicted_data_matrix : {np.isnan(predicted_data_matrix).sum()}")


    # Combine predicted data matrix with the existing data matrix, ensuring no duplication with 'ctrl'
    original_data_matrix = (
        original_data.X.toarray() if hasattr(original_data.X, "toarray") else original_data.X
    )
    
    # Step 4: Check for NaNs in the original data matrix
    print(f"Step 4: Nomber of NaN in the original_data_matrix : {np.isnan(original_data_matrix).sum()}")

    # Step 5: Combine matrices and check for NaNs in the combined data
    combined_data = np.vstack([original_data_matrix, predicted_data_matrix])
    print(f"Step 5: Nomber of NaN in combined_data  : {np.isnan(combined_data).sum()}")

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
    print(f"Step 6: Nomber of NaN in anndata_with_predictions.X : {np.isnan(anndata_with_predictions.X).sum()}")

    # Optional: If desired, set raw attribute for the AnnData object
    anndata_with_predictions.raw = anndata_with_predictions.copy()

    return anndata_with_predictions


# Load the dataset
dataset_path = "C:\\Users\\Shadow\\Desktop\\BioHack24\\scPRAM\\processed_datasets_all\\datasets\\scrna-lupuspatients\\kang-hvg.h5ad"
adata = sc.read_h5ad(dataset_path)
cell_types = adata.obs['cell_type'].unique()
output_dir = ".\\output_ood_models_1"

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
        'key' : 'cell_type',
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
    

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy import sparse
import numpy as np

# Initialisation du dictionnaire pour stocker les résultats R²
r2_results = {}

# Boucle sur chaque type cellulaire unique
for cell_type in anndata_with_predictions.obs['cell_type'].unique():
    # Filtrer les données pour ce type cellulaire
    cell_data = anndata_with_predictions[anndata_with_predictions.obs['cell_type'] == cell_type]

    # Extraire les données pour les conditions 'stim' et 'predicted'
    stim_data = cell_data[cell_data.obs['condition'] == 'stim'].X
    predicted_data = cell_data[cell_data.obs['condition'] == 'predicted'].X

    # Convertir les matrices creuses en matrices denses si nécessaire
    if sparse.issparse(stim_data):
        stim_data = stim_data.toarray()
    if sparse.issparse(predicted_data):
        predicted_data = predicted_data.toarray()

    # Calculer la moyenne d'expression pour chaque gène
    stim_mean = stim_data.mean(axis=0)
    predicted_mean = predicted_data.mean(axis=0)

    # Calculer la corrélation de Pearson
    r = np.corrcoef(stim_mean, predicted_mean)[0, 1]
    r2 = r ** 2  # R² basé sur la corrélation de Pearson

    # Stocker le résultat
    r2_results[cell_type] = r2
    print(f"R² pour {cell_type} entre les moyennes des gènes 'stim' et 'predicted' : {r2:.4f}")

    # Préparer les données pour la visualisation
    df_plot = pd.DataFrame({
        'Stim Mean Expression': stim_mean,
        'Predicted Mean Expression': predicted_mean
    })

    # Tracer le nuage de points avec la ligne de régression
    plt.figure(figsize=(8, 6))
    sns.regplot(
        x='Stim Mean Expression',
        y='Predicted Mean Expression',
        data=df_plot,
        scatter_kws={'s': 10},  # Taille des points
        line_kws={'color': 'red'}  # Couleur de la ligne de régression
    )
    plt.title(f'Regression Plot for {cell_type}\nR² = {r2:.4f}')
    plt.xlabel('Stim Mean Expression')
    plt.ylabel('Predicted Mean Expression')
    plt.grid(True)
    plt.show()

# Afficher tous les résultats R²
print("\nR² entre les moyennes d'expression pour chaque type cellulaire entre 'stim' et 'predicted':")
for cell_type, r2 in r2_results.items():
    print(f"{cell_type}: {r2:.4f}")
    


from scipy.spatial.distance import cdist
from scipy.sparse import issparse
import numpy as np
import pandas as pd
import scanpy as sc

def compute_edistance(set1, set2):
    """
    Compute the energy distance between two datasets.
    """
    intra_dist1 = np.mean(cdist(set1, set1, metric="euclidean"))
    intra_dist2 = np.mean(cdist(set2, set2, metric="euclidean"))
    inter_dist = np.mean(cdist(set1, set2, metric="euclidean"))
    return 2 * inter_dist - intra_dist1 - intra_dist2

def compute_perturbation_score_per_cell_type(anndata, 
                                             n_comps=50, 
                                             condition_col="condition", 
                                             stim_key="stim", 
                                             ctrl_key="ctrl", 
                                             pred_key="predicted"):
    """
    Compute the perturbation score for each cell type.

    Parameters:
        anndata: AnnData object containing gene expression data.
        n_comps: Number of principal components to use.
        condition_col: Column name in `.obs` that specifies the condition.
        stim_key: Key for the stimulated condition.
        ctrl_key: Key for the control condition.
        pred_key: Key for the predicted condition.
        
    Returns:
        A dictionary mapping each cell type to its perturbation score.
    """
    perturbation_scores = {}

    # Perform PCA on the data
    if n_comps > min(anndata.shape):
        n_comps = min(anndata.shape) - 1

    sc.tl.pca(anndata, svd_solver="arpack", n_comps=n_comps)
    print(f"PCA with {n_comps} components computed.\n")

    # Iterate over each cell type
    for cell_type in anndata.obs['cell_type'].unique():
        print(f"Processing cell type: {cell_type}")

        # Subset the data for the current cell type
        cell_data = anndata[anndata.obs['cell_type'] == cell_type]

        # Extract the subsets for stimulated, control, and predicted data
        stim_adata = cell_data[cell_data.obs[condition_col] == stim_key]
        ctrl_adata = cell_data[cell_data.obs[condition_col] == ctrl_key]
        pred_adata = cell_data[cell_data.obs[condition_col] == pred_key]

        # Skip if any subset is empty
        if stim_adata.shape[0] == 0 or ctrl_adata.shape[0] == 0 or pred_adata.shape[0] == 0:
            print(f"Skipping {cell_type} due to insufficient data.\n")
            continue

        # Extract PCA embeddings
        stim_pca = stim_adata.obsm["X_pca"]
        ctrl_pca = ctrl_adata.obsm["X_pca"]
        pred_pca = pred_adata.obsm["X_pca"]

        # Convert sparse matrices to dense
        if issparse(stim_pca): stim_pca = stim_pca.toarray()
        if issparse(ctrl_pca): ctrl_pca = ctrl_pca.toarray()
        if issparse(pred_pca): pred_pca = pred_pca.toarray()

        # Compute energy distances
        edistance_stim_pred = compute_edistance(stim_pca, pred_pca)  # Perturbed vs Predicted
        edistance_ctrl_pred = compute_edistance(ctrl_pca, pred_pca)  # Control vs Predicted

        # Avoid division by zero
        if edistance_ctrl_pred == 0:
            perturbation_score = np.nan
        else:
            perturbation_score = edistance_stim_pred / edistance_ctrl_pred

        perturbation_scores[cell_type] = perturbation_score
        print(f"Perturbation score for {cell_type}: {perturbation_score}\n")

    return perturbation_scores


perturbation_scores = compute_perturbation_score_per_cell_type(
    anndata=anndata_with_predictions,
    n_comps=50,
    condition_col="condition",
    stim_key="stim",
    ctrl_key="ctrl",
    pred_key="predicted"
)

# Display the results
print("Scaled perturbation scores for all cell types:")
for cell_type, score in perturbation_scores.items():
    print(f"{cell_type}: {score:.4f}")


#-------------- mmd -------------------

from sklearn.metrics.pairwise import polynomial_kernel, rbf_kernel
import numpy as np
import pandas as pd
import scanpy as sc

def compute_mmd(set1, set2, kernel="linear", **kernel_kwargs):
    """
    Compute the Maximum Mean Discrepancy (MMD) between two datasets.
    
    Parameters:
        set1: np.ndarray
            First dataset (e.g., real perturbed data).
        set2: np.ndarray
            Second dataset (e.g., predicted data).
        kernel: str
            Type of kernel to use. Options are 'linear', 'rbf', and 'poly'.
        **kernel_kwargs:
            Additional arguments for the kernel function (e.g., gamma for RBF).
            
    Returns:
        float
            MMD score.
    """
    if kernel == "linear":
        XX = np.dot(set1, set1.T)
        YY = np.dot(set2, set2.T)
        XY = np.dot(set1, set2.T)
    elif kernel == "rbf":
        XX = rbf_kernel(set1, set1, **kernel_kwargs)
        YY = rbf_kernel(set2, set2, **kernel_kwargs)
        XY = rbf_kernel(set1, set2, **kernel_kwargs)
    elif kernel == "poly":
        XX = polynomial_kernel(set1, set1, **kernel_kwargs)
        YY = polynomial_kernel(set2, set2, **kernel_kwargs)
        XY = polynomial_kernel(set1, set2, **kernel_kwargs)
    else:
        raise ValueError(f"Unsupported kernel type: {kernel}")

    return XX.mean() + YY.mean() - 2 * XY.mean()

def compute_mmd_per_cell_type(anndata, 
                              n_comps=50, 
                              condition_col="condition", 
                              stim_key="stim", 
                              ctrl_key="ctrl", 
                              pred_key="predicted", 
                              kernel="linear", 
                              **kernel_kwargs):
    """
    Compute the MMD for each cell type.

    Parameters:
        anndata: AnnData
            Annotated data matrix containing the data.
        n_comps: int
            Number of PCA components to use.
        condition_col: str
            Column in `.obs` specifying the condition of the cells.
        stim_key: str
            Key for the stimulated condition in `.obs`.
        ctrl_key: str
            Key for the control condition in `.obs`.
        pred_key: str
            Key for the predicted condition in `.obs`.
        kernel: str
            Kernel type for MMD. Options: 'linear', 'rbf', 'poly'.
        **kernel_kwargs:
            Additional parameters for the kernel function.

    Returns:
        dict
            A dictionary mapping each cell type to its MMD perturbation score.
    """
    mmd_scores = {}

    # Perform PCA on the data
    if n_comps > min(anndata.shape):
        n_comps = min(anndata.shape) - 1

    sc.tl.pca(anndata, svd_solver="arpack", n_comps=n_comps)
    print(f"PCA with {n_comps} components computed.\n")

    # Iterate over each cell type
    for cell_type in anndata.obs['cell_type'].unique():
        print(f"Processing cell type: {cell_type}")

        # Subset the data for the current cell type
        cell_data = anndata[anndata.obs['cell_type'] == cell_type]

        # Extract subsets for stimulated, control, and predicted data
        stim_adata = cell_data[cell_data.obs[condition_col] == stim_key]
        ctrl_adata = cell_data[cell_data.obs[condition_col] == ctrl_key]
        pred_adata = cell_data[cell_data.obs[condition_col] == pred_key]

        # Extract PCA embeddings
        stim_pca = stim_adata.obsm["X_pca"]
        ctrl_pca = ctrl_adata.obsm["X_pca"]
        pred_pca = pred_adata.obsm["X_pca"]

        # Compute MMD scores
        mmd_stim_pred = compute_mmd(stim_pca, pred_pca, kernel=kernel, **kernel_kwargs)  # Stimulated vs Predicted
        mmd_ctrl_pred = compute_mmd(ctrl_pca, pred_pca, kernel=kernel, **kernel_kwargs)  # Control vs Predicted

        # Combine scores into a perturbation score
        mmd_score = mmd_stim_pred / mmd_ctrl_pred
        mmd_scores[cell_type] = mmd_score
        print(f"MMD perturbation score for {cell_type}: {mmd_score}\n")

    return mmd_scores

# Compute MMD scores for all cell types
mmd_scores = compute_mmd_per_cell_type(
    anndata=anndata_with_predictions,
    n_comps=50,
    condition_col="condition",
    stim_key="stim",
    ctrl_key="ctrl",
    pred_key="predicted",
    kernel="linear"#,  # Example: using RBF kernel
#    gamma=1.0  # Example parameter for the RBF kernel
)

# Display the results
print("MMD perturbation scores for all cell types:")
for cell_type, score in mmd_scores.items():
    print(f"{cell_type}: {score:.4f}")

#------------------ euclidean distances -------------------

from sklearn.metrics.pairwise import euclidean_distances
import numpy as np
import pandas as pd
import scanpy as sc

def compute_mean_euclidean_distance(set1, set2):
    """
    Compute the mean Euclidean distance between two datasets.

    Parameters:
        set1: np.ndarray
            First dataset (e.g., real perturbed data).
        set2: np.ndarray
            Second dataset (e.g., predicted data).
            
    Returns:
        float
            Mean Euclidean distance between set1 and set2.
    """
    pairwise_distances = euclidean_distances(set1, set2)
    return pairwise_distances.mean()

def compute_euclidean_distance_per_cell_type(anndata, 
                                             n_comps=50, 
                                             condition_col="condition", 
                                             stim_key="stim", 
                                             ctrl_key="ctrl", 
                                             pred_key="predicted"):
    """
    Compute the mean Euclidean distance for each cell type.

    Parameters:
        anndata: AnnData
            Annotated data matrix containing the data.
        n_comps: int
            Number of PCA components to use.
        condition_col: str
            Column in `.obs` specifying the condition of the cells.
        stim_key: str
            Key for the stimulated condition in `.obs`.
        ctrl_key: str
            Key for the control condition in `.obs`.
        pred_key: str
            Key for the predicted condition in `.obs`.

    Returns:
        dict
            A dictionary mapping each cell type to its Euclidean perturbation score.
    """
    euclidean_scores = {}

    # Perform PCA on the data
    if n_comps > min(anndata.shape):
        n_comps = min(anndata.shape) - 1

    sc.tl.pca(anndata, svd_solver="arpack", n_comps=n_comps)
    print(f"PCA with {n_comps} components computed.\n")

    # Iterate over each cell type
    for cell_type in anndata.obs['cell_type'].unique():
        print(f"Processing cell type: {cell_type}")

        # Subset the data for the current cell type
        cell_data = anndata[anndata.obs['cell_type'] == cell_type]

        # Extract subsets for stimulated, control, and predicted data
        stim_adata = cell_data[cell_data.obs[condition_col] == stim_key]
        ctrl_adata = cell_data[cell_data.obs[condition_col] == ctrl_key]
        pred_adata = cell_data[cell_data.obs[condition_col] == pred_key]

        # Extract PCA embeddings
        stim_pca = stim_adata.obsm["X_pca"]
        ctrl_pca = ctrl_adata.obsm["X_pca"]
        pred_pca = pred_adata.obsm["X_pca"]

        # Compute Euclidean distances
        euclidean_stim_pred = compute_mean_euclidean_distance(stim_pca, pred_pca)  # Stimulated vs Predicted
        euclidean_ctrl_pred = compute_mean_euclidean_distance(ctrl_pca, pred_pca)  # Control vs Predicted

        # Combine scores into a perturbation score
        euclidean_score = euclidean_stim_pred / euclidean_ctrl_pred
        euclidean_scores[cell_type] = euclidean_score
        print(f"Euclidean perturbation score for {cell_type}: {euclidean_score}\n")

    return euclidean_scores

# Compute Euclidean distance scores for all cell types
euclidean_scores = compute_euclidean_distance_per_cell_type(
    anndata=anndata_with_predictions,
    n_comps=50,
    condition_col="condition",
    stim_key="stim",
    ctrl_key="ctrl",
    pred_key="predicted"
)

# Display the results
print("Euclidean perturbation scores for all cell types:")
for cell_type, score in euclidean_scores.items():
    print(f"{cell_type}: {score:.4f}")


# Bar plot with each metrics 


import matplotlib.pyplot as plt
import numpy as np

# Assurez-vous que les listes suivantes contiennent les données réelles générées
cell_types = list(r2_results.keys())  # Les types cellulaires


r2_scores = [float(val) for val in r2_results.values()]  # Conversion si nécessaire
edistances = [float(val) for val in perturbation_scores.values()]  # Conversion si nécessaire
mmd_res = [float(val) for val in mmd_scores.values()]  # Conversion des ArrayView
euclidean_dist = [float(val) for val in euclidean_scores.values()]  # Conversion si nécessaire




# Configuration des sous-graphiques
fig, axes = plt.subplots(1, 4, figsize=(16, 8), sharey=True)

# Graphique pour R² Scores
axes[0].barh(cell_types, r2_scores, color='blue', edgecolor='black')
axes[0].set_title("R² Scores")
axes[0].set_xlabel("Valeur")
axes[0].invert_yaxis()  # Alignement des types cellulaires sur tous les graphiques

# Graphique pour Energy Distance
axes[1].barh(cell_types, edistances, color='green', edgecolor='black')
axes[1].set_title("Energy Distance")
axes[1].set_xlabel("Valeur")

# Graphique pour MMD Scores
axes[2].barh(cell_types, mmd_res, color='orange', edgecolor='black')
axes[2].set_title("MMD Scores")
axes[2].set_xlabel("Valeur")

# Graphique pour Euclidean Distance
axes[3].barh(cell_types, euclidean_dist, color='red', edgecolor='black')
axes[3].set_title("Euclidean Distance")
axes[3].set_xlabel("Valeur")

# Ajuster l'espacement entre les sous-graphiques
plt.tight_layout()

# Afficher le graphique
plt.show()











