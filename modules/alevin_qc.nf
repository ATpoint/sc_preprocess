process AlevinQC {

    tag "$sample_id"

    errorStrategy 'finish'

    label 'process_alevinqc'

    publishDir params.outdir, mode: params.publishmode    

    container params.container

    input:
    tuple val(sample_id), path(alevin)
    val(suffix)

    output:
    path("*.html")
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:
    """
    Rscript --vanilla $baseDir/bin/alevin_qc.R $sample_id $alevin $suffix

    echo ${task.process}: > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt
    
    echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
    echo 'alevinQC:' \$(Rscript -e "cat(as.character(packageVersion('alevinQC')))") >> versions.txt
    """ 

}

