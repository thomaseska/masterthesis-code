#!/usr/bin/env/nextflow
nextflow.enable.dsl=2

process CREATE_DATASETS {
    input:
    each sc_count_mat
    each metadata
    each folder_path
    each dataset_name
    each method
    each spillover
    each overrep

    shell:
    '''
    ~/masterthesis/scripts/nextflow/create_datasets.R '!{sc_count_mat}' '!{metadata}' '!{folder_path}' '!{dataset_name}' '!{method}' '!{spillover}' '!{overrep}'
    '''
    
}

workflow {
    CREATE_DATASETS(params.sc_count_mat,
                    params.metadata,
                    params.folder_path,
                    params.dataset_name,
                    params.method,
                    params.spillover,
                    params.overrep)
}