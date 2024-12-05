# Perturbation Analysis Using scGen on Single-Cell RNA-Seq Data

---

## Introduction to scGen

**scGen** is a deep learning model designed to predict cellular responses to perturbations at the single-cell level. By leveraging variational autoencoders (VAEs), scGen captures the intrinsic variability within cell populations and models how cells respond to different stimuli or conditions. This enables researchers to simulate how specific cell types might react to new perturbations without the need for extensive wet-lab experiments.

**Key aspects of scGen**:
- **Generative Modeling**: scGen uses VAEs to learn a low-dimensional latent space representing the cellular state.
- **Perturbation Vector**: The model learns a perturbation vector that represents the shift in gene expression caused by a stimulus.
- **Prediction of Unseen Conditions**: scGen can predict how cells would respond to perturbations not present in the training data.
- **Cell-Type Specificity**: The model accounts for different cell types, allowing for cell-type-specific predictions.

---

## Pipeline Development Overview

To effectively utilize scGen for perturbation analysis, we developed a comprehensive analysis pipeline encapsulated within the `scGenPerturbationAnalysis` class. This class streamlines the workflow, making it modular and reusable for different datasets and experimental conditions.

### Main Advances and Steps

1. **Data Loading and Initialization**:
    - Loaded the **Stephenson et al. 2021** COVID-19 dataset using `pertpy`.
    - Initialized the `scGenPerturbationAnalysis` class with the loaded data.
2. **Data Exploration**:
    - Explored the dataset to understand the distribution of diseases and time points.
    - Visualized UMAP embeddings colored by various metadata to inspect data quality and clustering.
3. **Data Filtering**:
    - Filtered the data to include only healthy control samples to focus on specific perturbations.
    - Verified the filtering by checking counts and distributions post-filtering.
4. **Normalization and Visualization**:
    - Normalized the data to a total count per cell and (eventually) applied log transformation.
    - Visualized the data to ensure that normalization did not adversely affect the data structure.
5. **Data Preprocessing**:
    - Selected highly variable genes to reduce dimensionality and focus on informative genes.
    - Removed specific conditions (e.g., certain time points) and subsampled the data to manage computational load.
6. **Training Set Preparation**:
    - Excluded specific cell types under certain conditions from the training set to simulate an unseen condition during prediction.
7. **Model Setup and Training**:
    - Set up the AnnData object with appropriate batch and labels keys required by scGen.
    - Trained the scGen model using the prepared training set, adjusting hyperparameters as needed.
8. **Latent Space Visualization**:
    - Obtained latent representations from the trained model.
    - Visualized the latent space using UMAP to assess how well the model captures the data structure.
9. **Prediction**:
    - Used the trained scGen model to predict the stimulated state of specific cell types.
    - Assigned predicted labels to facilitate evaluation.
10. **Evaluation of Predictions**:
    - Combined real and predicted data for direct comparison.
    - Performed dimensionality reduction (PCA and UMAP) to visualize the overlap between predicted and real stimulated cells.
11. **Differential Expression Analysis**:
    - Identified differentially expressed genes between control and stimulated states.
    - Focused on top genes to understand the biological relevance of predictions.
12. **Correlation and Distance Metrics**:
    - Plotted mean gene expression correlation between predicted and actual stimulated cells.
    - Computed distance metrics (e.g., Earth Mover's Distance) to quantitatively assess prediction accuracy.

---

## Caveats and Issues Encountered

### Issue with Infinite Values During Preprocessing

**Problem**:
    - When setting `highly_variable_genes=True` in the `preprocess_data` method, we encountered the following error:
```
ValueError: cannot specify integer `bins` when input data contains infinity
```

**Cause**:
- The error was due to infinite values in the data when calculating means and dispersions for selecting highly variable genes.
- Not applying log transformation after normalization resulted in extremely large values, leading to infinite dispersions.

**Investigation**:
- Checked for infinite and NaN values in the data after normalization.
- Realized that without log transformation, the data retained a skewed distribution with extreme values.

### Additional Issues

- **Data Integrity**: Presence of zero or near-zero counts causing computational issues during dispersion calculations.
- **Zero Variance Genes**: Genes with zero variance can affect statistical calculations and model training.

---

## Solutions and Recommendations

1. **Apply Log Transformation**:
    - **Solution**: Enable log transformation in the `normalize_and_visualize` method by setting `log_transform=True`.
    - **Rationale**: Log transformation stabilizes variance and mitigates the effect of extreme values, preventing infinite dispersions.
```python
analysis.normalize_and_visualize(
    normalize=True,
    log_transform=True,  # Enable log transformation
    target_sum=1e4,
    plot_cols=['disease']
)
```
2. **Data Inspection and Cleaning**:
    - **Check for Infinite and NaN Values**:

```python
import numpy as np

print("Infinite values in data:", np.isinf(analysis.data.X).any())
print("NaN values in data:", np.isnan(analysis.data.X).any())
```
    - **Handle Infinite/NaN Values**:
        - Replace infinite values with NaN and remove affected genes or cells.
        - Filter out genes with zero variance to maintain data integrity.
3. **Ensure Data Integrity Before HVG Selection**:
    - Perform data scaling and filtering to remove outliers.
    - Verify that the data does not contain extreme values that could cause computational issues.
4. **Adjust Preprocessing Parameters**:
    - If issues persist, adjust parameters like `n_top_genes` or use alternative methods for highly variable gene selection.

---

## Conclusion

By addressing the issues encountered during preprocessing, we successfully implemented the perturbation analysis pipeline using scGen. The key takeaways and steps include:

- **Importance of Proper Data Preprocessing**: Ensuring data integrity through normalization and log transformation is crucial.
- **Modular Pipeline Design**: Encapsulating the workflow within a class enhances reusability and organization.
- **Thorough Evaluation**: Combining visual and quantitative methods provides a comprehensive assessment of model performance.
- **Flexibility in Methods**: Adjusting parameters and methods based on data characteristics is essential for successful analysis.

This pipeline serves as a robust framework for conducting perturbation analysis on single-cell RNA-seq data, enabling the exploration of cellular responses to various stimuli.

---

## References

- **scGen Paper**: Lotfollahi et al., ["scGen predicts single-cell perturbation responses"](https://www.nature.com/articles/s41592-019-0494-8)
- **Stephenson et al. 2021 Dataset**: [Single-cell multi-omics analysis of the immune response in COVID-19](https://www.nature.com/articles/s41591-021-01329-2)