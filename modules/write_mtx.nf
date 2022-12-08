process WriteMtx {

    tag "$basename"

    errorStrategy 'finish'

    label 'process_mtx'

    publishDir = [
        path: params.outdir, 
        mode: params.publishmode,
        saveAs: {filename -> filename.contains("ncells.txt") ? null : filename}
    ]

    container params.container

    input:
    tuple val(basename), path(quant_rna), path(quant_fb)
    path(expanded_features)
    path(gene2type)
    path(translate_table)
    
    output:
    path("${basename}*.mtx.gz"), emit: mtx
    path("${basename}_barcodes.tsv.gz"), emit: barcodes
    path("${basename}_features.tsv.gz"), emit: features
    path("${basename}_ncells.txt"), emit: ncells
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:
    
    def do_translate = params.translate_barcodes==true ? '--do_translate' : ''

    """
    Rscript --vanilla $baseDir/bin/mtx.R \
        --sampleid $basename \
        --alevin_rna $quant_rna \
        --alevin_fb $quant_fb \
        --expanded_features $expanded_features \
        --gene2type $gene2type \
        --translate_table $translate_table \
        $do_translate

    echo ${task.process}:${basename} > command_lines.txt
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