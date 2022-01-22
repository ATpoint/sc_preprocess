process WriteNcells {

    cpus 1
    memory 100.MB
    time '5min'

    publishDir = [
        path: params.outdir, 
        mode: params.publishmode,
    ]

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