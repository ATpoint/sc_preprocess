// Extract exonic and intronic (=unspliced) sequences from reference files using basically the R script from
// the COMBINE lab: https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/

process ParseExonIntronTx {

    cpus    params.threads
    memory  params.mem

    publishDir params.outdir, mode: params.pubmode

    input:
    path(genome)            
    path(gtf)               

    output:
    path "*annotation.expanded.fa.gz",              emit: txtome            // 1
    path "*annotation.expanded.gtf.gz",             emit: gtf               // 2
    path "*annotation.expanded.features.tsv.gz",    emit: features          // 3
    path "*annotation.expanded.tx2gene.tsv",        emit: tgmap             // 4
    path "*annotation.mtRNA.txt",                   emit: mtgenes           // 5
    path "*annotation.rRNA.txt",                    emit: rrnagenes         // 6
    path "html/*.html",                             emit: html              // 7
    path "*annotation.gene2type.txt",               emit: gene2type         // 8

    script:
    """
    Rscript --vanilla -e 'rmarkdown::render(input=\"${baseDir}/bin/parse_introns.rmd\", 
                                            output_file=\"parse_introns.html\",
                                            params=list(genome=\"${genome}\", 
                                                        gtf=\"${gtf}\",
                                                        gene_name=\"${params.gene_name}\", 
                                                        gene_type=\"${params.gene_type}\", 
                                                        gene_id=\"${params.gene_id}\",
                                                        readlength=\"${params.readlength}\",
                                                        chrM=\"${params.chrM}\"),
                                                        knit_root_dir=getwd(), 
                                                        output_dir=paste0(getwd(), \"/html\")
                                            )'
    """

}