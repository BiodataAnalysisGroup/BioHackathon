#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process CELL_ORACLE {
    container "${params.container__cell_oracle}"
    label 'high_proc'

    input:
        tuple val(tool_name), val(task_name), val(dataset_name), path(dataset_path)

    output:
        tuple val(tool_name), val(task_name), val(dataset_name), path("tool_output")

    when:
        tool_name == params.tool_cell_oracle

    script:
        """
            echo "Running tool Cell Oracle on dataset ${dataset_name}"
            mkdir tool_output
        """
}