// Read data into R, create spliced-unspliced count tables, save as mtx.gz module:

process WriteMtx {

    cpus    params.threads
    memory  params.mem

    publishDir params.outdir, mode: params.pubmode

    tag "Write mtx.gz: $basename"

    input:
    tuple val(basename), path(quant)
    path(expanded_features)
    path(gene2type)

    output:
    path("*.mtx.gz"), emit: mtx
    path("*.tsv.gz"), emit: colrowannot
    path("html/*.html")
    path("*.md"), optional: true

    script:
    """
    Rscript --vanilla -e 'rmarkdown::render(input=\"${baseDir}/bin/mtx.Rmd\", 
                                            output_file=\"${basename}_mtx.html\",
                                            params=list(expanded_features=\"${expanded_features}\", 
                                                        basename = \"${basename}\",
                                                        alevin_quant=\"${quant}\", 
                                                        gene2type=\"${gene2type}\"),
                                                        knit_root_dir=getwd(), 
                                                        output_dir=paste0(getwd(), \"/html\")
                                            )'
    """ 

}