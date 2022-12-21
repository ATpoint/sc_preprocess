process Quant {

    tag "${meta.id}"

    label 'process_quant_full'

    errorStrategy 'finish'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

    input:
    tuple val(meta), path(r1), path(r2)
    path(idx)
    path(tgmap)
    path(rrnagenes)
    path(mtgenes)

    output:
    tuple val(meta), path("$meta.id"), emit: quants
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    // as in https://combine-lab.github.io/alevin-tutorial/2020/alevin-features/
    script:
    """
    salmon alevin --no-version-check \
        -i $idx --tgMap $tgmap --mrna $mtgenes --rrna $rrnagenes \
        -o ${meta.id} -l $meta.libtype -p $task.cpus --dumpFeatures \
        $params.r1_type $params.r2_type $params.args -1 $r1 -2 $r2

    echo ${task.process}:${meta.id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'salmon:' \$(salmon --version | cut -d " " -f2) > versions.txt        
    """

}

process QuantFB {

    tag "${meta.id}"

    errorStrategy 'finish'

    label 'process_quant_features'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

    input:
    tuple val(meta), path(r1), path(r2)
    path(idx)
    path(tgmap)

    output:
    tuple val(meta), path("${meta.id}_fb"), emit: quants
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    // as in https://combine-lab.github.io/alevin-tutorial/2020/alevin-features/
    script:
    """  
    salmon alevin --no-version-check \
        -i $idx --tgMap $tgmap \
        -o ${meta.id}_fb -l $meta.libtype -p $task.cpus --dumpFeatures --keepCBFraction 1.0 \
        $params.r1_type $params.r2_type $params.args -1 $r1 -2 $r2

    echo ${task.process}:${meta.id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'salmon:' \$(salmon --version | cut -d " " -f2) > versions.txt                
    """

}