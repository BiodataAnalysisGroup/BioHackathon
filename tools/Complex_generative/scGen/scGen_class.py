import pertpy as pt
import scanpy as sc
from scvi import REGISTRY_KEYS
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class scGenPerturbationAnalysis:
    """
    A class to perform perturbation analysis using scGen on single-cell RNA-seq data.

    This class encapsulates the full analysis pipeline from data loading and preprocessing
    to model training, prediction, and evaluation.

    Attributes:
    - data: The original AnnData object containing the dataset.
    - train_set: The training dataset after excluding specific cells.
    - scgen_model: The trained scGen model.
    - latent_adata: AnnData object containing latent representations.
    - pred: AnnData object containing predicted cells.
    - delta: The perturbation vector learned by the model.
    - eval_adata: AnnData object used for evaluating predictions.
    - diff_genes: List of differentially expressed genes.
    - r2_value: R^2 value from mean expression correlation plot.
    - r_value_diff: R value for differentially expressed genes.
    - perturbation_score: Perturbation score computed using a distance metric.
    """

    def __init__(self, data):
        """
        Initialize the scGenPerturbationAnalysis class with the input data.

        Parameters:
        - data: AnnData object containing single-cell RNA-seq data.
        """
        self.data = data.copy()
        self.train_set = None
        self.scgen_model = None
        self.latent_adata = None
        self.pred = None
        self.delta = None
        self.eval_adata = None
        self.diff_genes = None
        self.r2_value = None
        self.r_value_diff = None
        self.perturbation_score = None

    def explore_data(self, groupby_cols=None, plot_cols=None, value_count_cols=None):
        """
        Perform basic exploration of the data, including printing counts, value counts, and plotting UMAPs.

        Parameters:
        - groupby_cols: List of column names to group by for counts (default: None).
        - plot_cols: List of column names to use for UMAP coloring (default: None).
        - value_count_cols: List of column names to display value counts (default: None).
        """
        print("Exploring data...")

        # Group by specified columns to understand sample distribution
        if groupby_cols is not None:
            combined_counts = self.data.obs.groupby(groupby_cols).size()
            combined_counts_df = combined_counts.reset_index(name='counts')
            print(f"Counts by {groupby_cols}:")
            print(combined_counts_df)
            print()
        else:
            print("No groupby columns provided for counts.")
            print()

        # Display value counts for specified columns
        if value_count_cols is not None:
            for col in value_count_cols:
                if col in self.data.obs.columns:
                    print(f"Value counts for '{col}':")
                    print(self.data.obs[col].value_counts())
                    print()
                else:
                    print(f"Column '{col}' not found in data.obs.")
        else:
            print("No columns provided for value counts.")
            print()

        # Visualize UMAP embeddings colored by specified annotations
        if plot_cols is not None:
            for col in plot_cols:
                if col in self.data.obs.columns:
                    sc.pl.umap(self.data, color=col, show=True)
                else:
                    print(f"Column '{col}' not found in data.obs.")
        else:
            print("No columns provided for UMAP plotting.")
        print()

    def filter_data(self, filter_conditions=None, verify_counts=False, count_groupby_cols=None):
        """
        Filter the data based on user-defined conditions.

        Parameters:
        - filter_conditions: Dictionary where keys are column names and values are the desired values to keep.
                            Values can be single values or lists of values.
                            Example: {'disease': ['normal', 'covid19']}
        - verify_counts: Boolean indicating whether to verify the filtering by checking counts (default: False).
        - count_groupby_cols: List of column names to group by when checking counts (default: None).
                            If None and verify_counts is True, defaults to filter_conditions.keys() if filter_conditions is provided.
                            If both are None, counts will not be printed.
        """
        print("Filtering data based on user-defined conditions...")
        if filter_conditions is not None:
            query_parts = []
            for col, vals in filter_conditions.items():
                if isinstance(vals, list):
                    # If vals is a list, create an OR condition for all values
                    vals_formatted = ', '.join([f'"{val}"' for val in vals])
                    query_parts.append(f'({col} in [{vals_formatted}])')
                else:
                    # Single value condition
                    query_parts.append(f'({col} == "{vals}")')
            query = " & ".join(query_parts)
            self.data = self.data[self.data.obs.eval(query)].copy()
            print(f"Filtered data with conditions: {filter_conditions}")
            print()
        else:
            print("No filter conditions provided. Data remains unfiltered.")
            print()

        # Verify the filtering by checking counts
        if verify_counts:
            if count_groupby_cols is None:
                if filter_conditions is not None:
                    count_groupby_cols = list(filter_conditions.keys())
                else:
                    print("No groupby columns provided for counts.")
                    count_groupby_cols = []
            if count_groupby_cols:
                print("Counts after filtering:")
                print(self.data.obs.groupby(count_groupby_cols).size().reset_index(name='counts'))
                print()
            else:
                print("No columns provided for grouping counts.")
                print()

    def normalize_and_visualize(self, normalize=True, log_transform=True, target_sum=1e4, plot_cols=None):
        """
        Normalize the data and visualize UMAP embeddings.

        Parameters:
        - normalize: Whether to normalize the data (default: True).
        - log_transform: Whether to log-transform the data (default: True).
        - target_sum: Target sum for normalization (default: 1e4).
        - plot_cols: List of column names to use for UMAP coloring (default: None).
        """
        if normalize:
            print("Normalizing data...")
            sc.pp.normalize_total(self.data, target_sum=target_sum)
            print(f"Data normalized to total counts per cell with target sum {target_sum}.")
            print()
        else:
            print("Normalization skipped.")
            print()
            
        if log_transform:
            sc.pp.log1p(self.data)
            print("Data log-transformed (log1p).")
            print()
        else:
            print("Log transformation skipped.")
            print()
            
        # Visualize UMAP embedding colored by specified columns
        if plot_cols is not None:
            for col in plot_cols:
                if col in self.data.obs.columns:
                    sc.pl.umap(self.data, color=col, show=True)
                else:
                    print(f"Column '{col}' not found in data.obs.")
        else:
            print("No columns provided for UMAP plotting.")
        print()

    def preprocess_data(
        self,
        highly_variable_genes=True,
        n_top_genes=2000,
        remove_conditions=None,
        subsample_n=None
    ):
        """
        Preprocess the data for scGen, including selecting highly variable genes,
        removing specific conditions, and subsampling.

        Parameters:
        - highly_variable_genes: Whether to select highly variable genes (default: True).
        - n_top_genes: Number of top variable genes to select (default: 2000).
        - remove_conditions: Dictionary of conditions to remove.
                             Example: {'time_after_LPS': ['90m']}
        - subsample_n: Number of observations to subsample (default: None).
        """
        print("Preprocessing data for scGen...")

        if highly_variable_genes:
            sc.pp.highly_variable_genes(self.data, n_top_genes=n_top_genes)
            print(f"Selected top {n_top_genes} highly variable genes.")
            print()
        else:
            print("Selection of highly variable genes skipped.")
            print()

        if remove_conditions is not None:
            for col, values in remove_conditions.items():
                if col in self.data.obs.columns:
                    self.data = self.data[~self.data.obs[col].isin(values)].copy()
                    print(f"Removed cells where '{col}' is in {values}.")
                else:
                    print(f"Column '{col}' not found in data.obs.")
            print()
        else:
            print("No conditions provided for removal.")
            print()

        if subsample_n is not None:
            sc.pp.subsample(self.data, n_obs=subsample_n)
            print(f"Subsampled data to {subsample_n} observations.")
            print()
        else:
            print("Subsampling skipped.")
            print()

        # Keep a copy of the raw data
        self.data.raw = self.data.copy()
        print("Preprocessing completed.")
        print()

    def prepare_training_set(self, exclude_query=None):
        """
        Prepare the training set by excluding specific cells based on a custom query.

        Parameters:
        - exclude_query: A string representing the query to exclude cells.
                        Example: '(cell_type == "B_naive") & (time_after_LPS == "10h")'
        """
        print("Preparing training set by excluding specific conditions...")
        if exclude_query is not None:
            self.train_set = self.data[~self.data.obs.eval(exclude_query)].copy()
            print(f"Excluded cells where: {exclude_query}")
            print()
        else:
            self.train_set = self.data.copy()
            print("No exclusion query provided. Using full dataset as training set.")
            print()

    def setup_anndata(self, batch_key=None, labels_key=None):
        """
        Set up the AnnData object for scGen.

        Parameters:
        - batch_key: Column name in obs to use as batch key (default: None).
        - labels_key: Column name in obs to use as labels key (default: None).
        """
        print("Setting up AnnData for scGen...")
        if batch_key is None or labels_key is None:
            print("Both 'batch_key' and 'labels_key' must be provided.")
            return
        pt.tl.Scgen.setup_anndata(
            self.train_set,
            batch_key=batch_key,
            labels_key=labels_key
        )
        print(f"AnnData set up for scGen with batch_key='{batch_key}' and labels_key='{labels_key}'.")
        print()

    def train_model(
        self,
        max_epochs=20,
        batch_size=32,
        early_stopping=True,
        accelerator="cpu",
        devices="auto",
        **kwargs
    ):
        """
        Initialize and train the scGen model.

        Parameters:
        - max_epochs: Maximum number of epochs to train (default: 20).
        - batch_size: Batch size for training (default: 32).
        - early_stopping: Whether to use early stopping (default: True).
        - accelerator: Device to use for training ('cpu' or 'gpu', default: 'cpu').
        - devices: Devices to use ('auto' or specific device, default: 'auto').
        - **kwargs: Additional keyword arguments for scGen training.
        """
        print("Initializing and training scGen model...")
        # Initialize the scGen model with the training data
        self.scgen_model = pt.tl.Scgen(self.train_set)
        print("scGen model initialized.")
        print()

        # Train the scGen model with specified parameters
        self.scgen_model.train(
            max_epochs=max_epochs,
            batch_size=batch_size,
            early_stopping=early_stopping,
            accelerator=accelerator,
            devices=devices,
            **kwargs
        )
        print("scGen model trained.")
        print()

        # Save the trained model
        self.scgen_model.save("model_perturbation_prediction.pt", overwrite=True)
        print("scGen model saved to 'model_perturbation_prediction.pt'.")
        print()

    def visualize_latent_space(self, plot_cols=None):
        """
        Get the latent representation from the model and visualize using UMAP.

        Parameters:
        - plot_cols: List of column names to use for UMAP coloring (default: None).
        """
        print("Visualizing latent space...")
        # Obtain the latent representations from the trained model
        latent_X = self.scgen_model.get_latent_representation()
        # Create a new AnnData object with the latent representations
        self.latent_adata = sc.AnnData(X=latent_X, obs=self.train_set.obs.copy())
        print("Latent representations obtained.")
        print()

        # Compute the neighborhood graph and UMAP embedding
        sc.pp.neighbors(self.latent_adata)
        sc.tl.umap(self.latent_adata)
        print("UMAP embedding computed.")
        print()

        # Plot the UMAP embedding colored by specified columns
        if plot_cols is not None:
            for col in plot_cols:
                if col in self.latent_adata.obs.columns:
                    sc.pl.umap(
                        self.latent_adata,
                        color=col,
                        wspace=0.4,
                        frameon=False,
                        show=True
                    )
                else:
                    print(f"Column '{col}' not found in latent_adata.obs.")
        else:
            print("No columns provided for UMAP plotting.")
        print("Latent space visualized.")
        print()

    def make_prediction(
        self,
        ctrl_key=None,
        stim_key=None,
        celltype_to_predict=None,
        condition_col=None,
        cell_type_col=None
    ):
        """
        Predict the stimulated state of specific cell types.

        Parameters:
        - ctrl_key: The control condition key (e.g., 'nan').
        - stim_key: The stimulated condition key (e.g., '10h').
        - celltype_to_predict: The cell type to predict (e.g., 'B_naive').
        - condition_col: The column name for conditions in obs (default: None).
        - cell_type_col: The column name for cell types in obs (default: None).
        """
        print("Making predictions...")
        if None in [ctrl_key, stim_key, celltype_to_predict, condition_col, cell_type_col]:
            print("All parameters must be provided for prediction.")
            return

        # Predict the stimulated state using the trained scGen model
        self.pred, self.delta = self.scgen_model.predict(
            ctrl_key=ctrl_key,
            stim_key=stim_key,
            celltype_to_predict=celltype_to_predict
        )
        # Assign a label to the predicted cells
        self.pred.obs[condition_col] = 'predicted'
        print("Predictions made.")
        print()

    def evaluate_prediction(
        self,
        condition_col=None,
        cell_type_col=None,
        control_condition=None,
        stimulated_condition=None,
        celltype_to_evaluate=None,
        pca_components=50
    ):
        """
        Evaluate the prediction by combining control, real stimulated, and predicted data.

        Parameters:
        - condition_col: The column name for conditions in obs (default: None).
        - cell_type_col: The column name for cell types in obs (default: None).
        - control_condition: The control condition key (e.g., 'nan').
        - stimulated_condition: The stimulated condition key (e.g., '10h').
        - celltype_to_evaluate: The cell type to evaluate (e.g., 'B_naive').
        - pca_components: Number of PCA components to compute (default: 50).
        """
        print("Evaluating predictions...")
        if None in [condition_col, cell_type_col, control_condition, stimulated_condition, celltype_to_evaluate]:
            print("All parameters must be provided for evaluation.")
            return

        # Extract control cells
        ctrl_adata = self.data[
            (self.data.obs[cell_type_col] == celltype_to_evaluate) &
            (self.data.obs[condition_col] == control_condition)
        ].copy()

        # Extract real stimulated cells
        stim_adata = self.data[
            (self.data.obs[cell_type_col] == celltype_to_evaluate) &
            (self.data.obs[condition_col] == stimulated_condition)
        ].copy()

        # Concatenate control, stimulated, and predicted data
        self.eval_adata = ctrl_adata.concatenate(stim_adata, self.pred)
        print("Evaluation data prepared.")
        print()

        # Perform PCA on the evaluation data
        sc.tl.pca(self.eval_adata, n_comps=pca_components)
        # Plot PCA colored by condition
        sc.pl.pca(
            self.eval_adata,
            color=condition_col,
            frameon=False,
            show=True
        )
        print("PCA plot created.")
        print()

        # Compute neighbors and UMAP for the evaluation data
        sc.pp.neighbors(self.eval_adata)
        sc.tl.umap(self.eval_adata)
        # Plot UMAP colored by condition
        sc.pl.umap(
            self.eval_adata,
            color=condition_col,
            frameon=False,
            show=True
        )
        print("UMAP plot created.")
        print()

    def identify_diff_genes(
        self,
        cell_type_col=None,
        condition_col=None,
        celltype_of_interest=None,
        method='wilcoxon'
    ):
        """
        Identify differentially expressed genes between control and stimulated cells.

        Parameters:
        - cell_type_col: The column name for cell types in obs (default: None).
        - condition_col: The column name for conditions in obs (default: None).
        - celltype_of_interest: The cell type to analyze (e.g., 'B_naive').
        - method: The method for differential expression (default: 'wilcoxon').
        """
        print("Identifying differentially expressed genes...")
        if None in [cell_type_col, condition_col, celltype_of_interest]:
            print("All parameters must be provided for identifying differentially expressed genes.")
            return

        # Select cells of the specified cell type
        cells_of_interest = self.data[self.data.obs[cell_type_col] == celltype_of_interest].copy()
        # Visualize UMAP of the cells
        sc.pl.umap(
            cells_of_interest,
            color=condition_col,
            frameon=False,
            show=True
        )
        print(f"UMAP of '{celltype_of_interest}' cells plotted.")
        print()

        # Perform differential expression analysis
        sc.tl.rank_genes_groups(cells_of_interest, groupby=condition_col, method=method)
        # Extract the names of differentially expressed genes
        groups = cells_of_interest.obs[condition_col].unique().tolist()
        if len(groups) >= 2:
            self.diff_genes = cells_of_interest.uns["rank_genes_groups"]["names"][groups[1]]
            print("Differentially expressed genes identified:")
            print(self.diff_genes)
        else:
            print("Not enough groups for differential expression analysis.")
        print()

    def plot_mean_correlation(
        self,
        condition_col=None,
        stimulated_condition=None,
        x_label='Predicted',
        y_label='Ground Truth',
        top_genes=10
    ):
        """
        Plot mean gene expression correlation between predicted and real stimulated cells.

        Parameters:
        - condition_col: The column name for conditions in obs (default: None).
        - stimulated_condition: The stimulated condition key (e.g., '10h').
        - x_label: Label for the x-axis (default: 'Predicted').
        - y_label: Label for the y-axis (default: 'Ground Truth').
        - top_genes: Number of top differentially expressed genes to highlight (default: 10).
        """
        print("Plotting mean gene expression correlation...")
        if None in [condition_col, stimulated_condition]:
            print("Both 'condition_col' and 'stimulated_condition' must be provided.")
            return

        # Get the original condition key used in scGen
        condition_key = self.scgen_model.adata_manager.get_state_registry(REGISTRY_KEYS.BATCH_KEY).original_key

        # Plot mean expression correlation for all genes
        self.r2_value = self.scgen_model.plot_reg_mean_plot(
            self.eval_adata,
            condition_key=condition_key,
            axis_keys={"x": "predicted", "y": stimulated_condition},
            labels={"x": x_label, "y": y_label},
            path_to_save="./reg_mean.pdf",
            show=True,
            legend=False
        )
        print(f"Mean expression correlation (R^2 value): {self.r2_value}")
        print()

        # Plot mean expression correlation for top differentially expressed genes
        if self.diff_genes is not None:
            self.r2_value, self.r_value_diff = self.scgen_model.plot_reg_mean_plot(
                self.eval_adata,
                condition_key=condition_key,
                axis_keys={"x": "predicted", "y": stimulated_condition},
                gene_list=self.diff_genes[:top_genes],
                top_100_genes=self.diff_genes,
                x_coeff=0.4,
                y_coeff=0.75,
                labels={"x": x_label, "y": y_label},
                path_to_save="./reg_mean_diff_genes.pdf",
                show=True,
                legend=False
            )
            print(f"Mean expression correlation for top {top_genes} genes (R^2 value): {self.r2_value}")
            print()
        else:
            print("No differentially expressed genes available for plotting.")
            print()

    def compute_distance_metric(
        self,
        condition_col=None,
        control_condition=None,
        stimulated_condition=None,
        pca_components=50,
        metric='edistance'
    ):
        """
        Compute a distance metric to evaluate prediction accuracy.

        Parameters:
        - condition_col: The column name for conditions in obs (default: None).
        - control_condition: The control condition key (e.g., 'nan').
        - stimulated_condition: The stimulated condition key (e.g., '10h').
        - pca_components: Number of PCA components to compute (default: 50).
        - metric: The distance metric to use ('edistance' or 'mmd', default: 'edistance').
        """
        print("Computing distance metric...")
        if None in [condition_col, control_condition, stimulated_condition]:
            print("All parameters must be provided for computing distance metric.")
            return

        # Perform PCA on the evaluation data
        sc.tl.pca(self.eval_adata, n_comps=pca_components)
        print(f"PCA with {pca_components} components computed.")
        print()

        # Extract the subsets for real stimulated, control, and predicted data
        stim_adata2 = self.eval_adata[self.eval_adata.obs[condition_col] == stimulated_condition]
        ctrl_adata2 = self.eval_adata[self.eval_adata.obs[condition_col] == control_condition]
        pred_adata2 = self.eval_adata[self.eval_adata.obs[condition_col] == 'predicted']

        # Initialize the distance metric
        distance_metric = pt.tools.Distance(metric=metric)
        print(f"Distance metric '{metric}' initialized.")
        print()

        # Extract PCA embeddings for each subset
        pert_pca = np.array(stim_adata2.obsm["X_pca"])  # Real perturbed data
        pred_pca = np.array(pred_adata2.obsm["X_pca"])  # Simulated perturbed data
        ctrl_pca = np.array(ctrl_adata2.obsm["X_pca"])  # Control data

        # Compute the perturbation score using the distance metric
        self.perturbation_score = distance_metric.compare_distance(
            pert=pert_pca,
            pred=pred_pca,
            ctrl=ctrl_pca,
            mode="simple"
        )
        print(f"Perturbation score computed: {self.perturbation_score}")
        print()

    def run_all(self, params):
        """
        Run all steps of the analysis in sequence using provided parameters.

        Parameters:
        - params: Dictionary containing parameters for each method.
        """
        self.explore_data(
            groupby_cols=params.get('explore_data', {}).get('groupby_cols'),
            plot_cols=params.get('explore_data', {}).get('plot_cols')
        )
        self.filter_data(
            filter_conditions=params.get('filter_data', {}).get('filter_conditions')
        )
        self.normalize_and_visualize(
            normalize=params.get('normalize_and_visualize', {}).get('normalize', True),
            target_sum=params.get('normalize_and_visualize', {}).get('target_sum', 1e4),
            plot_cols=params.get('normalize_and_visualize', {}).get('plot_cols')
        )
        self.preprocess_data(
            highly_variable_genes=params.get('preprocess_data', {}).get('highly_variable_genes', True),
            n_top_genes=params.get('preprocess_data', {}).get('n_top_genes', 2000),
            remove_conditions=params.get('preprocess_data', {}).get('remove_conditions'),
            subsample_n=params.get('preprocess_data', {}).get('subsample_n')
        )
        self.prepare_training_set(
            exclude_conditions=params.get('prepare_training_set', {}).get('exclude_conditions')
        )
        self.setup_anndata(
            batch_key=params.get('setup_anndata', {}).get('batch_key'),
            labels_key=params.get('setup_anndata', {}).get('labels_key')
        )
        self.train_model(
            **params.get('train_model', {})
        )
        self.visualize_latent_space(
            plot_cols=params.get('visualize_latent_space', {}).get('plot_cols')
        )
        self.make_prediction(
            ctrl_key=params.get('make_prediction', {}).get('ctrl_key'),
            stim_key=params.get('make_prediction', {}).get('stim_key'),
            celltype_to_predict=params.get('make_prediction', {}).get('celltype_to_predict'),
            condition_col=params.get('make_prediction', {}).get('condition_col'),
            cell_type_col=params.get('make_prediction', {}).get('cell_type_col')
        )
        self.evaluate_prediction(
            condition_col=params.get('evaluate_prediction', {}).get('condition_col'),
            cell_type_col=params.get('evaluate_prediction', {}).get('cell_type_col'),
            control_condition=params.get('evaluate_prediction', {}).get('control_condition'),
            stimulated_condition=params.get('evaluate_prediction', {}).get('stimulated_condition'),
            celltype_to_evaluate=params.get('evaluate_prediction', {}).get('celltype_to_evaluate'),
            pca_components=params.get('evaluate_prediction', {}).get('pca_components', 50)
        )
        self.identify_diff_genes(
            cell_type_col=params.get('identify_diff_genes', {}).get('cell_type_col'),
            condition_col=params.get('identify_diff_genes', {}).get('condition_col'),
            celltype_of_interest=params.get('identify_diff_genes', {}).get('celltype_of_interest'),
            method=params.get('identify_diff_genes', {}).get('method', 'wilcoxon')
        )
        self.plot_mean_correlation(
            condition_col=params.get('plot_mean_correlation', {}).get('condition_col'),
            stimulated_condition=params.get('plot_mean_correlation', {}).get('stimulated_condition'),
            x_label=params.get('plot_mean_correlation', {}).get('x_label', 'Predicted'),
            y_label=params.get('plot_mean_correlation', {}).get('y_label', 'Ground Truth'),
            top_genes=params.get('plot_mean_correlation', {}).get('top_genes', 10)
        )
        self.compute_distance_metric(
            condition_col=params.get('compute_distance_metric', {}).get('condition_col'),
            control_condition=params.get('compute_distance_metric', {}).get('control_condition'),
            stimulated_condition=params.get('compute_distance_metric', {}).get('stimulated_condition'),
            pca_components=params.get('compute_distance_metric', {}).get('pca_components', 50),
            metric=params.get('compute_distance_metric', {}).get('metric', 'edistance')
        )
        print("All analysis steps completed.")