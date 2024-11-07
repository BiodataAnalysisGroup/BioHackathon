import pertpy as pp
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

    def explore_data(self):
        """
        Perform basic exploration of the data, including printing counts and plotting UMAPs.
        """
        print("Exploring data...")
        # Group by 'disease' and 'time_after_LPS' to understand sample distribution
        combined_counts = self.data.obs.groupby(['disease', 'time_after_LPS']).size()
        combined_counts_df = combined_counts.reset_index(name='counts')
        print("Counts by disease and time_after_LPS:")
        print(combined_counts_df)
        print()

        # Group by 'disease', 'Status', and 'time_after_LPS' for more detailed counts
        combined_counts = self.data.obs.groupby(['disease', 'Status', 'time_after_LPS']).size()
        combined_counts_df = combined_counts.reset_index(name='counts')
        print("Counts by disease, Status, and time_after_LPS:")
        print(combined_counts_df)
        print()

        # Visualize UMAP embeddings colored by different annotations
        sc.pl.umap(self.data, color="disease", show=True)
        sc.pl.umap(self.data, color="time_after_LPS", show=True)
        sc.pl.umap(self.data, color="author_cell_type", show=True)

    def filter_data(self):
        """
        Filter the data to include only healthy controls (disease == 'normal').
        """
        print("Filtering data to include only healthy controls...")
        self.data = self.data[self.data.obs["disease"] == "normal"].copy()
        print("Filtered data to include only normal samples.")
        print()

        # Verify the filtering by checking counts
        combined_counts = self.data.obs.groupby(['disease', 'time_after_LPS']).size()
        combined_counts_df = combined_counts.reset_index(name='counts')
        print("Counts after filtering:")
        print(combined_counts_df)
        print()

    def visualize_data(self):
        """
        Normalize the data and visualize UMAP embeddings.
        """
        print("Normalizing data and visualizing UMAP embeddings...")
        # Normalize the data to total counts per cell
        sc.pp.normalize_total(self.data, target_sum=1e4)
        print("Data normalized to total counts per cell.")
        print()

        # Visualize UMAP embedding colored by 'disease'
        sc.pl.umap(self.data, color="disease", show=True)

        # Print counts of cells by 'time_after_LPS' and 'cell_type'
        combined_counts = self.data.obs.groupby(['time_after_LPS', 'cell_type']).size()
        combined_counts_df = combined_counts.reset_index(name='counts')
        print("Counts by time_after_LPS and cell_type:")
        print(combined_counts_df)
        print()

        # Display the counts of each cell type
        cell_type_counts = self.data.obs["cell_type"].value_counts()
        print("Cell type counts:")
        print(cell_type_counts)
        print()

    def preprocess_data(self,
                        apply_log1p = True):
        """
        Preprocess the data for scGen, including selecting highly variable genes,
        removing specific time points, and subsampling.
        """
        print("Preprocessing data for scGen...")
        # Identify the top 2000 highly variable genes
        #sc.pp.highly_variable_genes(self.data, n_top_genes=2000)
        #print("Selected top 2000 highly variable genes.")
        #print()

        # Remove cells at 90 minutes after LPS stimulation
        #self.data = self.data[self.data.obs["time_after_LPS"] != "90m"].copy()
        #print("Removed cells at 90 minutes after LPS stimulation.")
        #print()

        # Subsample the data to 3000 observations to reduce computational load
        #sc.pp.subsample(self.data, n_obs=3000)
        #print("Subsampled data to 3000 observations.")
        #print()

        # Display the counts of 'time_after_LPS' and 'author_cell_type' after subsampling
        #time_counts = self.data.obs["time_after_LPS"].value_counts()
        #cell_type_counts = self.data.obs["author_cell_type"].value_counts()
        #print("Time after LPS counts after subsampling:")
        #print(time_counts)
        #print()
        #print("Author cell type counts after subsampling:")
        #print(cell_type_counts)
        #print()

        # Keep a copy of the raw data
        self.data.raw = self.data.copy()
        # Normalize the data    
        sc.pp.normalize_total(self.data)
        print("Data normalized for scGen.")

        if apply_log1p == True:
            sc.pp.log1p(self.data)
            print("Log(x+1) applied to the data set for scGen.")
        
        print()

    def prepare_training_set(self,
                             condition_col = "time_after_LPS",
                             stim_key = "10h",
                             celltype_col = "author_cell_type",
                             celltype_to_predict = "B_naive"):
        """
        Remove cell type at stim from the training set to test generalization.
        """
        # Define the condition to exclude naive B cells at 10h
        condition = (
            (self.data.obs[condition_col] == stim_key) &
            (self.data.obs[celltype_col] == celltype_to_predict)
        )
        # Exclude the specified cells
        self.train_set = self.data[~condition].copy()
        print("Removed",  celltype_to_predict, "& ",stim_key," from the training set")
        print()

        # Display the counts of each cell type in the training set
        #cell_type_counts = self.train_set.obs["author_cell_type"].value_counts()
        #print("Author cell type counts in training set:")
        #print(cell_type_counts)
        print()

    def setup_anndata(self,condition_col,celltype_col):
        """
        Set up the AnnData object for scGen.
        """
        print("Setting up AnnData for scGen...")
        # Set up the training data with batch and label keys
        pp.tl.Scgen.setup_anndata(
            self.train_set,
            batch_key = condition_col,
            labels_key= celltype_col
        )
        print("AnnData set up for scGen.")
        print()

    def train_model(
        self,
        max_epochs=20,
        batch_size=32,
        early_stopping=True,
        accelerator="cpu",
        devices="auto"
    ):
        """
        Initialize and train the scGen model.

        Parameters:
        - max_epochs: Maximum number of epochs to train (default: 20).
        - batch_size: Batch size for training (default: 32).
        - early_stopping: Whether to use early stopping (default: True).
        - accelerator: Device to use for training ('cpu' or 'gpu', default: 'cpu').
        - devices: Devices to use ('auto' or specific device, default: 'auto').
        """
        print("Initializing and training scGen model...")
        # Initialize the scGen model with the training data
        self.scgen_model = pp.tl.Scgen(self.train_set)
        print("scGen model initialized.")
        print()

        # Train the scGen model with specified parameters
        self.scgen_model.train(
            max_epochs=max_epochs,
            batch_size=batch_size,
            early_stopping=early_stopping,
            accelerator=accelerator,
            devices=devices
        )
        print("scGen model trained.")
        print()

        # Save the trained model
        self.scgen_model.save("model_perturbation_prediction.pt", overwrite=True)
        print("scGen model saved to 'model_perturbation_prediction.pt'.")
        print()

    def visualize_latent_space(self):
        """
        Get the latent representation from the model and visualize using UMAP.
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

        # Plot the UMAP embedding colored by 'time_after_LPS' and 'author_cell_type'
        sc.pl.umap(
            self.latent_adata,
            color=["time_after_LPS", "author_cell_type"],
            wspace=0.4,
            frameon=False,
            show=True
        )
        print("Latent space visualized.")
        print()

    def make_prediction(self,
                        ctrl_key = "nan",
                        stim_key = "10h",
                        celltype_to_predict = "B_naive",
                        condition_col = "time_after_LPS"):
        """
        Predict the stimulated state of naive B cells at 10h.
        """
        print("Making predictions for ", celltype_to_predict, "& ", stim_key)
        # Predict the stimulated state using the trained scGen model
        self.pred, self.delta = self.scgen_model.predict(
            ctrl_key=ctrl_key,
            stim_key=stim_key,
            celltype_to_predict=celltype_to_predict
        )
        # Assign a label to the predicted cells
        self.pred.obs[condition_col] = "pred"
        print("Predictions made.")
        print()

    def evaluate_prediction(self,
                            celltype_col = "author_cell_type",
                            celltype_to_predict = "B_naive",
                            condition_col = "time_after_LPS",
                            ctrl_key = "nan", 
                            stim_key = "10h"
                            ):
        """
        Evaluate the prediction by combining control, real stimulated, and predicted data.
        """
        print("Evaluating predictions...")
        # Extract control naive B cells
        ctrl_adata = self.data[
            (self.data.obs[celltype_col] == celltype_to_predict) &
            (self.data.obs[condition_col] == ctrl_key)
        ].copy()

        # Extract real stimulated naive B cells at 10h
        stim_adata = self.data[
            (self.data.obs[celltype_col] == celltype_to_predict) &
            (self.data.obs[condition_col] == stim_key)
        ].copy()

        # Concatenate control, stimulated, and predicted data
        self.eval_adata = ctrl_adata.concatenate(stim_adata, self.pred)
        print("Evaluation data prepared.")
        print()

        # Perform PCA on the evaluation data
        sc.tl.pca(self.eval_adata)
        # Plot PCA colored by 'time_after_LPS'
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
        # Plot UMAP colored by 'time_after_LPS'
        sc.pl.umap(
            self.eval_adata,
            color=condition_col,
            frameon=False,
            show=True
        )
        print("UMAP plot created.")
        print()

    def identify_diff_genes(self,
                            celltype_col = "author_cell_type",
                            celltype_to_predict = "B_naive",
                            condition_col = "time_after_LPS",
                            stim_key = "10h"):

        """
        Identify differentially expressed genes between control and stimulated B_naive cells.
        """
        print("Identifying differentially expressed genes...")
        # Select naive B cells from the original data
        celltype = self.data[self.data.obs[celltype_col] ==celltype_to_predict].copy()
        # Perform differential expression analysis
        sc.tl.rank_genes_groups(celltype, groupby=condition_col, method="wilcoxon")
        # Extract the names of differentially expressed genes for '10h' vs 'nan'
        self.diff_genes = celltype.uns["rank_genes_groups"]["names"][stim_key]
        print("Differentially expressed genes identified:")
        print(self.diff_genes)
        print()

    def plot_mean_correlation(self,
                              stim_key = "10h"):
        """
        Plot mean gene expression correlation between predicted and real stimulated cells.
        """
        print("Plotting mean gene expression correlation...")
        # Get the original condition key used in scGen
        condition_key = self.scgen_model.adata_manager.get_state_registry(REGISTRY_KEYS.BATCH_KEY).original_key

        # Plot mean expression correlation for all genes
        self.r2_value = self.scgen_model.plot_reg_mean_plot(
            self.eval_adata,
            condition_key=condition_key,
            axis_keys={"x": "pred", "y": stim_key},
            labels={"x": "Predicted", "y": "Ground Truth"},
            path_to_save="./reg_mean.pdf",
            show=True,
            legend=False
        )
        print(f"Mean expression correlation (R^2 value): {self.r2_value}")
        print()

        # Plot mean expression correlation for top differentially expressed genes (this should be modified and R^2 calculated
        # only for diff expressed genes)
        #self.r2_value, self.r_value_diff = self.scgen_model.plot_reg_mean_plot(
        #    self.eval_adata,
        #    condition_key=condition_key,
        #    axis_keys={"x": "pred", "y": stim_key},
        #    gene_list=self.diff_genes[:10],
        #    top_100_genes=self.diff_genes,
        #    x_coeff=3,
        #    y_coeff=1,
        #    labels={"x": "Predicted", "y": "Ground Truth"},
        #    path_to_save="./reg_mean_diff_genes.pdf",
        #    show=True,
        #    legend=False
        #)
        #print(f"Mean expression correlation for top genes (R^2 value): {self.r2_value}")
        print()

    def compute_distance_metric(self,
                                n_comps = 50,
                                metric = "edistance",
                                condition_col = "time_after_LPS",
                                stim_key = "10h",
                                ctrl_key = "nan"
                                ):
        """
        Compute a distance metric (e.g., energy distance) to evaluate prediction accuracy.
        """
        print("Computing distance metric...")
        if n_comps > min(self.eval_adata.shape):
            n_comps=min(self.eval_adata.shape) - 1

        sc.tl.pca(self.eval_adata, n_comps)
        print("PCA with", n_comps, "components computed.")
        print()

        # Extract the subsets for real stimulated, control, and predicted data
        stim_adata2 = self.eval_adata[self.eval_adata.obs[condition_col] == stim_key]
        ctrl_adata2 = self.eval_adata[self.eval_adata.obs[condition_col] == ctrl_key]
        pred_adata2 = self.eval_adata[self.eval_adata.obs[condition_col] == "pred"]

        # Initialize the distance metric (e.g., energy distance)
        distance_metric = pp.tools.Distance(metric)
        print("Distance metric initialized.")
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

    def run_all(self):
        """
        Run all steps of the analysis in sequence.
        """
        self.explore_data()
        self.filter_data()
        self.visualize_data()
        self.preprocess_data()
        self.prepare_training_set()
        self.setup_anndata()
        self.train_model()
        self.visualize_latent_space()
        self.make_prediction()
        self.evaluate_prediction()
        self.identify_diff_genes()
        self.plot_mean_correlation()
        self.compute_distance_metric()
        print("All analysis steps completed.")