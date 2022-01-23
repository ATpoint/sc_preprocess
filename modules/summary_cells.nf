process WriteNcells {

    cpus 1
    memory 2.GB

    errorStrategy 'finish'

    publishDir = [
        path: params.outdir, 
        mode: params.publishmode,
    ]

    if(workflow.profile.contains('conda'))  { conda "$params.environment" }
    if(workflow.profile.contains('docker')) { container "$params.container" }
    if(workflow.profile.contains('singularity')) { container "$params.container" }

    input:
    path(dirs)

    output:
    path("summary_detected_cells.txt")
    path("summary_detected_cells.pdf")

    script:
    """
    Rscript --vanilla $baseDir/bin/summary.R
    """

}