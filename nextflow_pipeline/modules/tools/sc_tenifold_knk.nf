#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process SC_TENIFOLD_KNK {
    container "${params.container__sc_tenifold_knk}"
    label 'high_proc'

    input:
        tuple val(tool_name), val(task_name), val(dataset_name), path(dataset_path)

    output:
        tuple val(tool_name), val(task_name), val(dataset_name), val(perturbed_gene), path(perturbed_genes_path), path(count_matrix_path)

    when:
        tool_name == params.tool_sc_tenifold_knk

    script:
        """
            echo "Running tool scTenifoldKnk on dataset ${dataset_name}"
            echo "Creating count matrix CSV file..."
            python ${projectDir}/../tools/GRN-prioritizing/scTenifoldKnk/extractCountMatrix.py \
                --sc_data ${dataset_path} \
                --out_path ${count_matrix_path}
            echo "Creating KO perturbation..."
            Rscript ${projectDir}/../tools/GRN-prioritizing/scTenifoldKnk/scTenifoldKnk.py \
                --countMatrix ${count_matrix_path} \
                --perturbedGenes ${perturbed_genes_path} \
                --geneKnockout ${perturbed_gene}   
        """
}