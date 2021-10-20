// Read data into R, create spliced-unspliced count tables, save as mtx.gz module:

process WriteMtx {

    cpus    params.threads
    memory  params.mem

    publishDir params.outdir, mode: params.pubmode

    tag "mtx.gz: $basename"

    input:
    tuple val(basename), path(quant)
    path(expanded_features)
    path(gene2type)

    output:
    path("*.mtx.gz"), emit: mtx
    path("*.tsv.gz"), emit: colrowannot
    
    script:
    """
    Rscript --vanilla $baseDir/bin/mtx.R \"${basename}\" ${quant}  ${expanded_features} ${gene2type} 
    """ 

}