// Read data into R, create spliced-unspliced count tables, save as mtx.gz module:

process WriteMtx {

    tag "$basename"

    errorStrategy 'finish'

    label 'process_mtx'

    publishDir = [
        path: params.outdir, 
        mode: params.publishmode,
        saveAs: {filename -> filename.contains("ncells.txt") ? null : filename}
    ]

    if(workflow.profile.contains('conda'))  { conda "$params.environment" }
    if(workflow.profile.contains('docker')) { container "$params.container" }
    if(workflow.profile.contains('singularity')) { container "$params.container" }

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
    
    script:
    // The quant_
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
    """ 

}