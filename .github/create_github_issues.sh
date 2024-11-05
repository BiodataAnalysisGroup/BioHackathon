#!/bin/bash

# Path to guidelines
GUIDELINES=$(cat /Users/marinaesteban-medina/Desktop/BioHackathon/.github/git_guidelines.txt)

# Create issues for dataset preprocessing
# gh issue create --title "Preprocess Dataset for Generative AI Models" \
# --body "$GUIDELINES\n\n### Task: Preprocess dataset for Generative AI models\n- **Description:** Create a preprocessing script to convert dataset IDs into a common format, including metadata needed for each perturbation task and metric calculation.\n- **Input:** Dataset ID\n- **Output:** Common formatted dataset with metadata.\n\n**Tags:** #python, #nextflow" \
# --label "python" --label "nextflow"

# gh issue create --title "Preprocess Dataset for GRN Prioritizing Models" \
# --body "$GUIDELINES\n\n### Task: Preprocess dataset for GRN prioritizing models\n- **Description:** Create a preprocessing script to format datasets with necessary metadata for GRN analysis and compatibility with downstream metrics.\n- **Input:** Dataset ID\n- **Output:** Common formatted dataset with metadata.\n\n**Tags:** #python, #nextflow" \
# --label "python" --label "nextflow"

# Create issues for Generative AI tools
generative_tools=("scGen https://github.com/theislab/scgen"
                  "scPRAM https://github.com/jiang-q19/scPRAM"
                  "scVIDR https://github.com/BhattacharyaLab/scVIDR"
                  "cellOT https://github.com/bunnech/cellot"
                  "trVAE https://github.com/theislab/trVAE")

for tool_info in "${generative_tools[@]}"; do
    tool_name=$(echo "$tool_info" | cut -d ' ' -f 1)
    tool_url=$(echo "$tool_info" | cut -d ' ' -f 2)
    
    gh issue create --title "Implement Tool Script for $tool_name" \
    --body "$(printf "$GUIDELINES\n\n### Task: Implement script for [$tool_name]($tool_url)\n- **Description:** Develop a script or set of scripts to run $tool_name. Identify required inputs and expected outputs. These scripts should ensure compatibility with Nextflow modules.\n- **Input:** Dataset data, metadata, and specified perturbation task\n- **Output:** Prediction results from $tool_name\n\n**Tags:** #python, #GenerativeAI")" \
    --label "python" --label "GenerativeAI" --label "tool"
done

# Create issues for GRN tools
grn_tools=("genKI https://github.com/yjgeno/GenKI"
           "scTenifold https://github.com/cailab-tamu/scTenifoldKnk"
           "D-Spin https://github.com/JialongJiang/DSPIN")

for tool_info in "${grn_tools[@]}"; do
    tool_name=$(echo "$tool_info" | cut -d ' ' -f 1)
    tool_url=$(echo "$tool_info" | cut -d ' ' -f 2)
    
    gh issue create --title "Implement Tool Script for $tool_name" \
    --body "$(printf "$GUIDELINES\n\n### Task: Implement script for [$tool_name]($tool_url)\n- **Description:** Create a script or set of scripts to run $tool_name, defining input data and output format. This should prepare the tool's output for integration into Nextflow.\n- **Input:** Dataset data, metadata, and perturbation task\n- **Output:** Prediction results from $tool_name\n\n**Tags:** #python, #GRNs")" \
    --label "python" --label "GRN-prioritizing" --label "tool"
done

# Create issues for metrics
# gh issue create --title "Calculate Distribution Distance Metrics" \
# --body "$GUIDELINES\n\n### Task: Calculate set of distribution distance metrics using pertpy\n- **Description:** Create scripts to calculate distribution distance metrics. Ensure compatibility with outputs from tools and necessary inputs for downstream metrics.\n- **Reference:** [Pertpy distribution metrics](https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/distances.html)\n- **Location:** /Users/marinaesteban-medina/Desktop/BioHackathon/.github\n\n**Tags:** #python, #metrics, #nextflow" \
# --label "python" --label "metrics" --label "nextflow"

# gh issue create --title "Calculate R2 and Pearson Correlation Metrics" \
# --body "$GUIDELINES\n\n### Task: Calculate R2 metric and Pearson correlation\n- **Description:** Develop scripts to calculate R2 and Pearson correlation metrics. Ensure they are compatible with tool outputs and meet metric input requirements.\n\n**Tags:** #python, #metrics, #nextflow" \
# --label "python" --label "metrics" --label "nextflow"

# gh issue create --title "Calculate AUC for GRN Reconstruction" \
# --body "$GUIDELINES\n\n### Task: Calculate AUC metric for GRN reconstruction\n- **Description:** Develop scripts to compute AUC for GRN reconstruction tasks, tailored to process tool outputs and required metric inputs.\n\n**Tags:** #python, #metrics, #nextflow" \
# --label "python" --label "metrics" --label "nextflow"

# # Create issues for Nextflow modules
# gh issue create --title "Nextflow Module for Dataset Preprocessing" \
# --body "$GUIDELINES\n\n### Task: Create Nextflow module for dataset preprocessing\n- **Description:** Develop a Nextflow module that utilizes the preprocessing script for dataset formatting and metadata integration across perturbation tasks.\n\n**Tags:** #nextflow" \
# --label "nextflow"

# gh issue create --title "Nextflow Module for Tool Scripts" \
# --body "$GUIDELINES\n\n### Task: Create Nextflow module for tool scripts\n- **Description:** Design a Nextflow module to standardize input/output for all tool scripts, preparing data for downstream analysis.\n\n**Tags:** #nextflow" \
# --label "nextflow"

# gh issue create --title "Nextflow Module for Metric Calculation" \
# --body "$GUIDELINES\n\n### Task: Create Nextflow module for metric calculations\n- **Description:** Develop a Nextflow module to run metric calculations, integrating with outputs from tool scripts and formatting for result interpretation.\n\n**Tags:** #nextflow" \
# --label "nextflow"
