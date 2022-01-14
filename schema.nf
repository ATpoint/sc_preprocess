#! /usr/bin/env nextflow

/* 
        SCHEMA DEFINITION FOR PARAMS VALIDATION FOR SCPREPROCESS
*/

def Map schema = [:] // don't change this line

// --------------------------------------------------------------------------------------------------------------

overall_outdir        = "$launchDir/sc_preprocess_results/"

// generic options:
schema.min_nf_version = [value: '21.10.6', type: 'string', mandatory: true, allowed: '']
schema.publishmode    = [value: 'copy', type: 'string', mandatory: true, allowed:['symlink', 'rellink', 'link', 'copy', 'copyNoFollow', 'move']]
schema.outdir         = [value: overall_outdir, type: 'string', mandatory: true]

// samplesheet:
schema.samplesheet    = [value: '', type: 'string', pattern: /.*\.csv$/]

// related to indexing the gentrome
schema.genome         = [value: 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz',
                         type: 'string', mandatory: true]
schema.gtf            = [value: 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz',
                         type: 'string', mandatory: true]
schema.gene_name      = [value: 'gene_name', type: 'string', mandatory: true]                 
schema.gene_id        = [value: 'gene_id', type: 'string', mandatory: true]
schema.gene_type      = [value: 'gene_type', type: 'string', mandatory: true]
schema.chrM           = [value: 'chrM', type: 'string', mandatory: true]
schema.rrna           = [value: 'rRNA', type: 'string', mandatory: true]
schema.idx_outdir     = [value: "${overall_outdir}/alevinIndex/", type: 'string', mandatory: true]
schema.idx_name       = [value: 'idx_gentrome', type: 'string', mandatory: true]
schema.idx_args       = [value: '--sparse --gencode', type: 'string']

// related to indexing the feature barcodes
schema.features_file  = [value: 'i', type: 'string', mandatory: false]
schema.features_name  = [value: 'idx_features', type: 'string', mandatory: false]

// related to quantification of reads against gentrome or feature barcode library
schema.quant_outdir   = [value: "${overall_outdir}/alevinQuant/", type: 'string']
schema.quant_libtype  = [value: 'ISR', type: 'string']
schema.quant_args     = [value: '--chromiumV3', type: 'string']
schema.quant_sf_args  = [value: '--bc-geometry 1[1-16] --umi-geometry 1[17-28] --read-geometry 2[11-25] --keepCBFraction 1.0', type: 'string'] // chromiumV3 with totalSeqB (15nt in R2)
schema.quant_sf_sfx   = [value: '_SF', type: 'string']

// related to mtx
schema.mtx_outdir     = [value: "${overall_outdir}/mtx/", type: 'string', mandatory: true]    

// related to the container/environment for the R/Bioconductor part of this workflow
schema.container      = [value:'atpoint/sc_preprocess:v1.3.2', type:'string', mandatory:true]
schema.environment    = [value: "$baseDir/environment.yml", type:'string', mandatory: true ]

// --------------------------------------------------------------------------------------------------------------

return schema // don't change this line