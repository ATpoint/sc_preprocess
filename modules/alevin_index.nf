
process AlevinIndex {

    label 'process_idx'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "bioconda::salmon=1.6.0"}
    if(workflow.profile.contains('docker')) { container "quay.io/biocontainers/salmon:1.6.0--h84f40af_0" }
    if(workflow.profile.contains('singularity')) { container "quay.io/biocontainers/salmon:1.6.0--h84f40af_0" }

    input:
    path(genome) 
    path(txtome) 
    val(idxname)
        
    output:
    path("decoynames.txt")
    path("gentrome.fa.gz")
    path(idxname), emit: idx
    
    script: 

    def decoynames  = "decoynames.txt"
    def gentrome    = "gentrome.fa.gz"

    """
    gzip -cd $genome | grep '^>' | cut -d " " -f 1 | awk '{gsub(\">\",\"\");print}' > $decoynames

    cat $txtome $genome > $gentrome

    salmon index --no-version-check \
        -t $gentrome \
        -d $decoynames \
        -i $idxname \
        -p $task.cpus \
        $params.additional
    """                

}

process AlevinIndexSF {

    cpus 1
    memory 100.MB

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "bioconda::salmon=1.6.0"}
    if(workflow.profile.contains('docker')) { container "quay.io/biocontainers/salmon:1.6.0--h84f40af_0" }
    if(workflow.profile.contains('singularity')) { container "quay.io/biocontainers/salmon:1.6.0--h84f40af_0" }

    input:
    path(featureset) 
    val(idxname)
        
    output:
    path(idxname), emit: idx
    path("$idxname/tgmap.txt"), emit: tgmap
    
    script: 
    """
    salmon index --no-version-check \
        -t ${featureset} \
        -i ${idxname} \
        --features \
        -k7

    awk 'OFS=\"\t\" {print \$1, \$1}' $featureset > "$idxname/tgmap.txt"
    """                

}