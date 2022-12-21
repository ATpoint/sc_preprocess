process AlevinQC {

    tag "${meta.id}"

    errorStrategy 'finish'

    label 'process_alevinqc'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]    

    container params.container

    input:
    tuple val(meta), path(quants)

    output:
    path("*.html")
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:
    
    """
    Rscript --vanilla $baseDir/bin/alevin_qc.R ${meta.id} $quants

    echo ${task.process}:${meta.id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt
    
    echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
    echo 'alevinQC:' \$(Rscript -e "cat(as.character(packageVersion('alevinQC')))") >> versions.txt
    """ 

}