process WriteNcells {

    cpus 1
    memory 2.GB

    errorStrategy 'finish'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]
    container params.container

    input:
    path(dirs)

    output:
    path("summary_detected_cells.txt")
    path("summary_detected_cells.pdf")
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions

    script:
    """
    Rscript --vanilla $baseDir/bin/summary.R

    echo ${task.process}: > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt
    
    echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
    """

}