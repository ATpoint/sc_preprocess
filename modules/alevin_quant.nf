// Alevin quant module:

process AlevinQuant {

    tag "$sample_id"

    label params.label

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(reads) // 1
    path(idx)                         // 2
    path(tgmap)                       // 3
    path(rrnagenes)                   // 4
    path(mtgenes)                     // 5

    output:
    path("${sample_id}"), emit: quants
    
    script:
    """
    salmon alevin --no-version-check \
        -i $idx --tgMap $tgmap --mrna $mtgenes --rrna $rrnagenes \
        -o $sample_id \
        -l $params.libtype \
        -p $task.cpus \
        $params.additional \
        -1 ${reads[0]} -2 ${reads[1]}
    """

}