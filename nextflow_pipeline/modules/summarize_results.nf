#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process SUMMARIZE_RESULTS {
    container "${params.container__summarize}"
    label 'low_proc'

    input:
        tuple val(tool_name_array), val(task_name), val(database_name_array), val(metric_name_array), path(metric_output_array)

    output:
        path("summary")

    script:
        """           
            mkdir summary
        """

}