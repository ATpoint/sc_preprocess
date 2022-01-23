// Read data into R, create spliced-unspliced count tables, save as mtx.gz module:

process AlevinQC {

    tag "$sample_id"

    errorStrategy 'finish'

    label 'process_alevinqc'

    publishDir params.outdir, mode: params.publishmode    

    if(workflow.profile.contains('conda'))  { conda "$params.environment" }
    if(workflow.profile.contains('docker')) { container "$params.container" }
    if(workflow.profile.contains('singularity')) { container "$params.container" }

    input:
    tuple val(sample_id), path(alevin)
    val(suffix)

    output:
    path("*.html")
    
    script:
    """
    Rscript --vanilla $baseDir/bin/alevin_qc.R $sample_id $alevin $suffix
    """ 

}

