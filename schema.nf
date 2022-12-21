#! /usr/bin/env nextflow

/*
 *
 *  SCHEMA DEFINITION FOR PARAMS VALIDATION
 *
 */

def Map schema = [:] // don't change this line

// --------------------------------------------------------------------------------------------------------------

schema.title1         = [title: 'GENERAL OPTIONS']
schema.min_nf_version = [value: '22.10.4', type: 'string', mandatory: true, allowed: '']
schema.publishmode    = [value: 'copy', type: 'string', mandatory: true, allowed:['symlink', 'rellink', 'link', 'copy', 'copyNoFollow', 'move']]
schema.outdir         = [value: "$launchDir/sc_preprocess_results/", type: 'string', mandatory: true]
schema.idx_outdir     = [value: "${schema.outdir['value']}/alevin_idx/", type: 'string', mandatory: true]
schema.quant_outdir   = [value: "${schema.outdir['value']}/alevin_quant/", type: 'string']
schema.mtx_outdir     = [value: "${schema.outdir['value']}/mtx/", type: 'string', mandatory: true]
schema.pipe_dir       = [value: "${schema.outdir['value']}/pipeline_info/", type: 'string', mandatory: true]
schema.qc_outdir      = [value: "${schema.outdir['value']}/qc_dir/", type: 'string', mandatory: true]
schema.container      = [value:'atpoint/sc_preprocess:v1.7.0', type:'string', mandatory:true]
schema.environment    = [value: "$baseDir/environment.yml", type:'string', mandatory: true ]
schema.samplesheet    = [value: '', type: 'string', pattern: /.*\.csv$/]

schema.title2         = [title: 'INDEXING OPTIONS']
schema.only_idx       = [value: false, type: 'logical']
schema.idx            = [value: '', type: 'string', mandatory: false]
schema.genome         = [value: '', type: 'string', mandatory: false]
schema.gtf            = [value: '', type: 'string', mandatory: false]
schema.gene_id        = [value: 'gene_id', type: 'string', mandatory: true]
schema.gene_name      = [value: 'gene_name', type: 'string', mandatory: true]                 
schema.read_length    = [value: 150, type: 'numeric', mandatory: true]
schema.gene_type      = [value: 'gene_type', type: 'string', mandatory: true]
schema.chrM           = [value: 'chrM', type: 'string', mandatory: true]
schema.rrna           = [value: 'rRNA', type: 'string', mandatory: true]
schema.idx_name       = [value: 'idx_gentrome', type: 'string',]
schema.idx_args       = [value: '--sparse', type: 'string']
schema.idx            = [value: '', type: 'string']
schema.tgmap          = [value: '', type: 'string']
schema.rrnagenes      = [value: '', type: 'string']
schema.mtgenes        = [value: '', type: 'string']
schema.features       = [value: '', type: 'string']
schema.gene2type      = [value: '', type: 'string']

schema.title3             = [title: 'QUANTIFICATION OPTIONS']
schema.r1_type            = [value: '--bc-geometry 1[1-16] --umi-geometry 1[17-28]', type: 'string', mandatory: true]
schema.r2_type            = [value: '--read-geometry 2[1-91]', type: 'string', mandatory: true]
schema.r2_type_fb         = [value: '--read-geometry 2[11-25]', type: 'string', mandatory: true]
schema.quant_args         = [value: '', type: 'string']
schema.quant_args_fb      = [value: '', type: 'string']
schema.features_file      = [value: '', type: 'string', mandatory: false]
schema.translate_barcodes = [value: false, type: 'logical']
schema.translate_list     = [value: "", type: 'string']

// --------------------------------------------------------------------------------------------------------------

return schema // don't change this line
