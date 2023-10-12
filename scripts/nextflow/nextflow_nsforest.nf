#!/usr/bin/env/nextflow
nextflow.enable.dsl=2

process RUN_NSFOREST {
    shell:
    '''
    python3 ~/masterthesis/scripts/ns_forest_signature.py 
    '''
}

workflow {
    RUN_NSFOREST()
}