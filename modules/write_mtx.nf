process Mtx {

    tag "$meta.id"

    errorStrategy 'finish'

    label 'process_mtx'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || 
                              filename.equals("command_lines.txt") ||
                              filename.equals("ncells.txt") ? null : filename } 
    ]

    container params.container

    input:
    tuple val(meta), path(quants)
    path(expanded_features)
    path(gene2type)
    path(translate_table)
    
    output:
    path("${meta.id}*.mtx.gz"), emit: mtx
    path("${meta.id}_barcodes.tsv.gz"), emit: barcodes
    path("${meta.id}_features.tsv.gz"), emit: features
    path("${meta.id}_ncells.txt"), emit: ncells
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:
    def with_fb = quants instanceof List // If a feature barcode quant is present quants is a List, else a TaskPath
    if(!with_fb){

        """
        Rscript --vanilla $baseDir/bin/mtx.R --sampleid ${meta.id} --alevin_rna $quants --expanded_features $expanded_features --gene2type $gene2type

        echo ${task.process}:${meta.id} > command_lines.txt
        cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt
        
        echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
        echo 'fishpond:' \$(Rscript -e "cat(as.character(packageVersion('fishpond')))") >> versions.txt
        echo 'Matrix:' \$(Rscript -e "cat(as.character(packageVersion('Matrix')))") >> versions.txt
        echo 'optparse:' \$(Rscript -e "cat(as.character(packageVersion('optparse')))") >> versions.txt
        echo 'rtracklayer:' \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))") >> versions.txt
        echo 'tximeta:' \$(Rscript -e "cat(as.character(packageVersion('tximeta')))") >> versions.txt
        echo 'SummarizedExperiment:' \$(Rscript -e "cat(as.character(packageVersion('SummarizedExperiment')))") >> versions.txt
        """ 

    } else {

        def do_translate = params.translate_barcodes ? "--translate_table $translate_table --do_translate" : ''

        """
        Rscript --vanilla $baseDir/bin/mtx.R \
            --sampleid ${meta.id} --alevin_rna ${quants[0]} --alevin_fb ${quants[1]} \
            --expanded_features $expanded_features --gene2type $gene2type $do_translate

        echo ${task.process}:${meta.id} > command_lines.txt
        cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt
        
        echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
        echo 'fishpond:' \$(Rscript -e "cat(as.character(packageVersion('fishpond')))") >> versions.txt
        echo 'Matrix:' \$(Rscript -e "cat(as.character(packageVersion('Matrix')))") >> versions.txt
        echo 'optparse:' \$(Rscript -e "cat(as.character(packageVersion('optparse')))") >> versions.txt
        echo 'rtracklayer:' \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))") >> versions.txt
        echo 'tximeta:' \$(Rscript -e "cat(as.character(packageVersion('tximeta')))") >> versions.txt
        echo 'SummarizedExperiment:' \$(Rscript -e "cat(as.character(packageVersion('SummarizedExperiment')))") >> versions.txt
        """ 

    }

}