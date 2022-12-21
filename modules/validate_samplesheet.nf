process ValidateSamplesheet {

    cpus 1
    memory 2.GB
    time '5m'

    errorStrategy 'finish'

     publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") || filename.equals("samplesheet_validated.csv") ? null : filename } 
    ]

    container params.container

    input:
    path(samplesheet) 
        
    output:
    path("samplesheet_validated.txt"), emit: samplesheet
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
        
    script: 
    """
    Rscript --vanilla $baseDir/bin/validate_samplesheet.R $samplesheet
    
    echo ${task.process}: > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt
    
    echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
    """                

}