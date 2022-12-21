process ParseExonIntronTx {

    // Implements part of: https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/
    
    label 'process_parse_exon_intron'

    errorStrategy 'finish'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

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
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions

    script:
    """
    Rscript --vanilla $baseDir/bin/parse_introns.R $genome $gtf $params.gene_id $params.gene_name $params.gene_type $params.read_length $params.chrM $params.rrna

    echo ${task.process}: > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt
    
    echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
    echo 'Biostrings:' \$(Rscript -e "cat(as.character(packageVersion('Biostrings')))") >> versions.txt
    echo 'BSgenome:' \$(Rscript -e "cat(as.character(packageVersion('BSgenome')))") >> versions.txt
    echo 'eisaR:' \$(Rscript -e "cat(as.character(packageVersion('eisaR')))") >> versions.txt
    echo 'GenomicFeatures:' \$(Rscript -e "cat(as.character(packageVersion('GenomicFeatures')))") >> versions.txt
    echo 'rtracklayer:' \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))") >> versions.txt
    """

}