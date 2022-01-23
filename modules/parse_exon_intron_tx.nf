// Extract exonic and intronic (=unspliced) sequences from reference files using basically the R script from
// the COMBINE lab: https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/

process ParseExonIntronTx {

    label 'process_parse_exon_intron'

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "$params.environment" }
    if(workflow.profile.contains('docker')) { container "$params.container" }
    if(workflow.profile.contains('singularity')) { container "$params.container" }

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
    path "*annotation.gene2type.txt",               emit: gene2type         // 7

    script:
    """
    Rscript --vanilla \
        $baseDir/bin/parse_introns.R \
            $genome $gtf \"$params.gene_name\" \"$params.gene_id\" \"$params.gene_type\" \"$params.chrM\" \"$params.rrna\"
    """

}