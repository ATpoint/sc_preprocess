process AlevinQuant {

    tag "$sample_id"

    label 'process_quant_full'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "bioconda::salmon=1.6.0"}
    if(workflow.profile.contains('docker')) { container "quay.io/biocontainers/salmon:1.6.0--h84f40af_0" }
    if(workflow.profile.contains('singularity')) { container "quay.io/biocontainers/salmon:1.6.0--h84f40af_0" }

    input:
    tuple val(sample_id), path(R1), path(R2), val(type)
    path(idx)                         // 2
    path(tgmap)                       // 3
    path(rrnagenes)                   // 4
    path(mtgenes)                     // 5

    output:
    tuple val(sample_id), path(sample_id), emit: quants
    
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
    """

}

process AlevinQuantFB {

    tag "$sample_id"

    label 'process_quant_features'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "bioconda::salmon=1.6.0"}
    if(workflow.profile.contains('docker')) { container "quay.io/biocontainers/salmon:1.6.0--h84f40af_0" }
    if(workflow.profile.contains('singularity')) { container "quay.io/biocontainers/salmon:1.6.0--h84f40af_0" }

    input:
    tuple val(sample_id), path(R1), path(R2), val(notused)
    path(idx)                         
    path(tgmap)     

    output:
    tuple val(sample_id), path("${sample_id}${params.suffix}"), emit: quants
    
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
        -1 ${R1} -2 ${R2}
    """

}