// Read data into R, create spliced-unspliced count tables, save as mtx.gz module:

process WriteMtx {

    tag "$basename"

    label 'process_mtx'

    publishDir params.outdir, mode: params.publishmode    

    if(workflow.profile.contains('conda'))  { conda "$params.environment" }
    if(workflow.profile.contains('docker')) { container "$params.container" }
    if(workflow.profile.contains('singularity')) { container "$params.container" }

    input:
    tuple val(basename), path(quant)
    path(expanded_features)
    path(gene2type)

    output:
    path("*.mtx.gz"), emit: mtx
    path("*barcodes.tsv.gz"), emit: barcodes
    path("*features.tsv.gz")
    
    script:
    """
    Rscript --vanilla $baseDir/bin/mtx.R \"${basename}\" ${quant}  ${expanded_features} ${gene2type} 
    """ 

}

process WriteMtxSF {

    tag "$basename"

    label 'process_mtx'

    publishDir params.outdir, mode: params.publishmode    

    if(workflow.profile.contains('conda'))  { conda "$params.environment" }
    if(workflow.profile.contains('docker')) { container "$params.container" }
    if(workflow.profile.contains('singularity')) { container "$params.container" }

    input:
    tuple val(basename), path(quant)

    output:
    path("*.mtx.gz"), emit: mtx
    
    script:
    """
    Rscript --vanilla $baseDir/bin/mtx_sf.R \"${basename}\" ${quant}
    """ 

}