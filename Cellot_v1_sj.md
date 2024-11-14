# Cellot Model Training and Evaluation: OOD Workflow

This document provides a step-by-step guide to preparing data, training the Cellot model in Out-of-Distribution (OOD) mode, and evaluating model predictions. This process involves custom modifications to the Cellot codebase to address specific requirements and improve model functionality (loss ouputs and anndata with ctrl/stim/pred).

---

## 1. Data Preparation

The Cellot model requires an AnnData object that contains information for two conditions:
  - **Control condition** (e.g., `ctrl`)
  - **Perturbed condition** (e.g., `stim`)

### Data Requirements
- **Data format**: The data should be in an AnnData structure compatible with single-cell analysis tools like Scanpy. Each observation in AnnData should include metadata in the `.obs` attribute, including cell type and condition.
  
- **Data Normalization**: 
  - Normalization is essential for consistent model performance. Using the `normalize_total` and `log1p` function in Scanpy.
  - **Scaling**: its specific impact on Cellot is still being evaluated, though standardizing features across cells may be beneficial for model performance..don't know yet

After preparing and normalizing data, save the AnnData object for input into the Cellot model training and evaluation scripts.

---

## 2. Training the Model with `cellot_train_v3_ood.py`

### Environment Setup
Follow the environment setup instructions from the Cellot GitHub repository to ensure dependencies are properly installed. Specifically, a Conda environment is recommended for managing dependencies.

### Custom Code Modifications
To handle specific issues encountered during model training, some modifications were made to the Cellot source code. These adjustments enhance compatibility with my OOD training and include:
  - **Files modified**: 
    - `cellot.data.cell`
    - `cellot.models.cellot`
    - `cellot.networks.icnns`
    - `cellot.train.train`

### Training Configuration
The Cellot OOD training script (`cellot_train_v3_ood.py`) includes a loop to automatically train individual models for each cell type. The key training parameters include:
- **Condition column** (`condition`): Defines the grouping of data into control and perturbed conditions.
- **Source and target conditions**: These specify the training setup. For example, `source='ctrl'` and `target='stim'`.
- **Epochs and batch size**: Standard parameters for deep learning models.
- **Holdout cell type** (`datasplit_holdout`): Specifies the cell type excluded from training for each OOD model.
  
### Running the Training Script
Run the `cellot_train_v3_ood.py` script after verifying all dependencies and data requirements. This script will:
1. Train the model for each specified cell type in OOD mode, using other cell types as training data.
2. Save models and training outputs, including loss tracking.

**Note**: Loss curves for transport functions are recorded and can be plotted at the end of each training session, though further integration into the loop is in progress.

---

## 3. Model Evaluation with `cellot_eval_v3_ood.py`

Once models are trained, the `cellot_eval_v3_ood.py` script enables evaluation of each cell-type-specific model. Evaluation includes visualizing predictions and computing performance metrics.

### Evaluation Outputs
1. **Dimensionality Reduction (PCA and UMAP)**:
   - The script generates PCA and UMAP visualizations, specifically for the holdout cell type (excluded during training) for each trained model.
   - These plots allow direct visual inspection of predicted cell distributions compared to actual data, providing insights into the model's performance in the OOD setting.

2. **Performance Metrics (in progress)**:
   - **R² Score**: Calculating R² for the predicted versus actual values helps quantify the model’s prediction accuracy for each cell type.
   - **Transport Distance**: In progress !! The distances (euclidian, e and mmd), a transport function metric, assesses how accurately the model translates control cells into their perturbed states.

The evaluation process allows detailed analysis of each model's performance per cell type, facilitating further adjustments and optimization of model parameters.

---

### Summary of Parameters and Model Configurations

- **`condition`**: Defines the control and target conditions in the dataset.
- **`datasplit_mode`**: When set to `ood`, this parameter splits the data to ensure the holdout cell type is excluded from the training data.
- **`datasplit_groupby`**: Controls grouping for the data split. For example, setting `['celltype','condition']` splits based on both cell type and condition.
  
---

This document provides the foundational steps for leveraging Cellot’s OOD training mode. Adapting the model to specific datasets and further optimizing parameters will enhance performance, especially with custom data configurations. It's still ugly with redundancy.
