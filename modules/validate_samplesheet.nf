process ValidateSamplesheet {

    cpus 1
    memory 2.GB
    time '5m'

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "$params.environment" }
    if(workflow.profile.contains('docker')) { container "$params.container" }
    if(workflow.profile.contains('singularity')) { container "$params.container" }

    input:
    path(samplesheet) 
    val(basedir)
    val(launchdir)
    val(projectdir)
        
    output:
    path("samplesheet_validated.txt"), emit: ssheet
        
    script: 
    """
    Rscript --vanilla $baseDir/bin/validate_samplesheet.R $samplesheet $basedir $launchdir $projectdir
    """                

}