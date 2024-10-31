#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process PREPROCESS_DATASET {
    container "${params.container__preprocess_dataset}"
    label 'low_proc'

    input:
        tuple val(dataset), path(input_dir)

    output:
        tuple val(dataset), path("dataset")

    script:
        """
            echo "Preprocessing dataset ${dataset.name} from ${input_dir}"
            mkdir dataset
        """

}