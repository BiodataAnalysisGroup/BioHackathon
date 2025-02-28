# scPRAM: Predicting Single-Cell Gene Expression Responses

scPRAM is a computational method designed to predict single-cell gene expression responses to various perturbations using attention mechanisms. It addresses challenges in obtaining perturbed samples and the high costs associated with large-scale single-cell sequencing experiments.

## Key Features

- **Alignment of Cell States**: Utilizes variational autoencoders and optimal transport to align cellular states before and after perturbation, facilitating accurate predictions.

- **Attention Mechanism**: Employs attention mechanisms to predict gene expression responses in unseen cell types by focusing on extracted information from the data.

## Performance Highlights

- **Ease of Implementation**: The tutorial is up-to-date and runs smoothly, making it easy to implement.

- **Interoperability**: Compatible with other tools like scGEN, enhancing its utility in various workflows.

- **Cross-Platform Availability**: Available for all operating systems; tested on Mac and Linux.

- **Sensitivity to Dataset Naming**: Functions are sensitive to the names of the predicted and training datasets; ensure consistency with tutorial examples.

## User Feedback

Users have reported that scPRAM provides accurate predictions and is easy to integrate into existing pipelines. The tool's compatibility with scGEN has been particularly appreciated, allowing for seamless integration into various workflows. However, some users have noted that the functions are sensitive to the names of the predicted and training datasets, requiring careful attention to naming conventions.

## Availability

The scPRAM tool is open-source and available at [https://github.com/jiang-q19/scPRAM](https://github.com/jiang-q19/scPRAM). For citation purposes, refer to the DOI: [https://doi.org/10.1093/bioinformatics/btae265](https://doi.org/10.1093/bioinformatics/btae265).
