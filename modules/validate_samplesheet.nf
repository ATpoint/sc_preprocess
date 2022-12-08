process ValidateSamplesheet {

    cpus 1
    memory 2.GB
    time '5m'

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    container params.container

    input:
    path(samplesheet) 
    val(basedir)
    val(launchdir)
    val(projectdir)
        
    output:
    path("samplesheet_validated.txt"), emit: ssheet
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
        
    script: 
    """
    Rscript --vanilla $baseDir/bin/validate_samplesheet.R $samplesheet $basedir $launchdir $projectdir

    echo ${task.process}: > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt
    
    echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
    """                

}