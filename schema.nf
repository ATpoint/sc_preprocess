#! /usr/bin/env nextflow

/*
 *
 *  SCHEMA DEFINITION FOR PARAMS VALIDATION
 *
 */

def Map schema = [:] // don't change this line

// --------------------------------------------------------------------------------------------------------------

// generic options:
schema.min_nf_version = [value: '21.10.6', type: 'string', mandatory: true, allowed: '']
schema.publishmode    = [value: 'copy', type: 'string', mandatory: true, allowed:['symlink', 'rellink', 'link', 'copy', 'copyNoFollow', 'move']]
schema.outdir         = [value: "$launchDir/sc_preprocess_results/", type: 'string', mandatory: true]
schema.pipe_dir       = [value: "${schema.outdir['value']}/pipeline_info/", type: 'string', mandatory: true]

// samplesheet:
schema.samplesheet    = [value: '', type: 'string', pattern: /.*\.csv$/]

// related to indexing the gentrome
schema.idx_only       = [value: false, type: 'logical']
schema.genome         = [value: 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz',
                         type: 'string', mandatory: true]
schema.gtf            = [value: 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz',
                         type: 'string', mandatory: true]
schema.gene_name      = [value: 'gene_name', type: 'string', mandatory: true]                 
schema.gene_id        = [value: 'gene_id', type: 'string', mandatory: true]
schema.gene_type      = [value: 'gene_type', type: 'string', mandatory: true]
schema.chrM           = [value: 'chrM', type: 'string', mandatory: true]
schema.rrna           = [value: 'rRNA', type: 'string', mandatory: true]
schema.idx_outdir     = [value: "${schema.outdir['value']}/alevin_idx/", type: 'string', mandatory: true]
schema.idx_name       = [value: 'idx_gentrome', type: 'string', mandatory: true]
schema.idx_args       = [value: '--sparse --gencode', type: 'string']

// related to providing an existing index:
schema.idx               = [value: '', type: 'string']
schema.tgmap             = [value: '', type: 'string']
schema.rrnagenes         = [value: '', type: 'string']
schema.mtrnagenes        = [value: '', type: 'string']
schema.expanded_features = [value: '', type: 'string']
schema.gene2type         = [value: '', type: 'string']

// related to indexing the feature barcodes
schema.features_file  = [value: '', type: 'string', mandatory: false]

// related to quantification of reads against gentrome or feature barcode library
schema.r1_type            = [value: '--bc-geometry 1[1-16] --umi-geometry 1[17-28]', type: 'string', mandatory: true]
schema.r2_type            = [value: '--read-geometry 2[1-91]', type: 'string', mandatory: true]
schema.r2_type_fb         = [value: '--read-geometry 2[11-25]', type: 'string', mandatory: true]
schema.quant_outdir       = [value: "${schema.outdir['value']}/alevin_quant/", type: 'string']
schema.libtype            = [value: 'ISR', type: 'string']
schema.quant_args         = [value: '', type: 'string']
schema.quant_args_fb      = [value: '', type: 'string']
schema.fb_suffix          = [value: '_fb', type: 'string']
schema.translate_barcodes = [value: true, type: 'logical', mandatory: true]
schema.translate_list     = [value: "$baseDir/assets/3M-february-2018.txt.gz", type: 'string', mandatory: true]

// related to mtx
schema.mtx_outdir     = [value: "${schema.outdir['value']}/mtx/", type: 'string', mandatory: true]
schema.rna_suffix     = [value: '_rna', type: 'string']

// related to alevinQC:
schema.qc_outdir      = [value: "${schema.outdir['value']}/alevin_qc/", type: 'string', mandatory: true]

// related to the container/environment for the R/Bioconductor part of this workflow
schema.container      = [value:'atpoint/sc_preprocess:v1.6.1', type:'string', mandatory:true]
schema.environment    = [value: "$baseDir/environment.yml", type:'string', mandatory: true ]

// --------------------------------------------------------------------------------------------------------------

return schema // don't change this line
