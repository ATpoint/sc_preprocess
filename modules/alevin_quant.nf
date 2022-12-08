process AlevinQuant {

    tag "$sample_id"

    errorStrategy 'finish'

    label 'process_quant_full'

    publishDir params.outdir, mode: params.publishmode

    container params.container

    input:
    tuple val(sample_id), path(R1), path(R2), val(type)
    path(idx)                         // 2
    path(tgmap)                       // 3
    path(rrnagenes)                   // 4
    path(mtgenes)                     // 5

    output:
    tuple val(sample_id), path(sample_id), emit: quants
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:
    """
    salmon alevin --no-version-check \
        -i $idx --tgMap $tgmap --mrna $mtgenes --rrna $rrnagenes \
        -o $sample_id \
        -l $params.libtype \
        -p $task.cpus \
        $params.additional \
        --dumpFeatures \
        -1 ${R1} -2 ${R2}

    echo ${task.process}:${sample_id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'salmon:' \$(salmon --version | cut -d " " -f2) > versions.txt        
    """

}

process AlevinQuantFB {

    tag "$sample_id"

    errorStrategy 'finish'

    label 'process_quant_features'

    publishDir params.outdir, mode: params.publishmode

    container params.container

    input:
    tuple val(sample_id), path(R1), path(R2), val(notused)
    path(idx)                         
    path(tgmap)     

    output:
    tuple val(sample_id), path("${sample_id}${params.suffix}"), emit: quants
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    // as in https://combine-lab.github.io/alevin-tutorial/2020/alevin-features/
    script:
    """   
    salmon alevin --no-version-check \
        -i $idx --tgMap $tgmap \
        -o ${sample_id}${params.suffix} \
        --libType $params.libtype \
        -p $task.cpus \
        $params.additional \
        --dumpFeatures \
        --keepCBFraction 1.0 \
        -1 ${R1} -2 ${R2}

    echo ${task.process}:${sample_id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'salmon:' \$(salmon --version | cut -d " " -f2) > versions.txt                
    """

}