#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process COSINE {
    container "${params.container__metrics}"
    label 'single_proc'

    input:
        tuple val(tool_name), val(task_name), val(dataset_name), val(metric_name), path(output_tool)

    output:
        tuple val(tool_name), val(task_name), val(dataset_name), val(metric_name), path("${metric_name}_${tool_name}_${dataset_name}.csv")

    when:
        metric_name == params.metric_cosine

    script:
        """
            echo "Calculating ${metric_name} for ${tool_name} predicting ${task_name} on dataset ${dataset_name}"
            touch ${metric_name}_${tool_name}_${dataset_name}.csv
        """
}