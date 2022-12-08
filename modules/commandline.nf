process CommandLines {

    cpus 1
    memory 500.MB

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    container params.container

    input:
    path(commands, stageAs: "?/*")
    path(versions, stageAs: "?/*")
    
    output:
    path("command_lines.txt")
    path("software_versions.txt")
        
    script:
    """
    find . -type d \
    | while read p; do
        echo \$(head -n 1 \${p}/command_lines.txt) >> commands.txt
        paste <(echo '') <(cat \${p}/command_lines.txt | sed '1d') >> commands.txt
      done
      
    Rscript --vanilla $baseDir/bin/sort_commands.R commands.txt       
    cat $versions | awk NF | sort --ignore-case -u > software_versions.txt
    """

}
