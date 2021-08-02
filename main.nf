#! /usr/bin/env nextflow

/*  
    ---------------------------------------------------------------------------------------------------------------
    Workflow for 10X scRNA-seq:

    --||    Creation of an index with salmon comprising the spliced- and unspliced transcriptome
            with the entire genome as decoy, starting from GENCODE reference files.

    --||    Quantification of the fastq files with Salmon/Alevin.

    --||    Read into R using tximeta, create spliced- and unspliced count tables, save as mtx.gz 
    
    ---------------------------------------------------------------------------------------------------------------
*/ 

nextflow.enable.dsl=2

def latest_sha = "git rev-parse HEAD".execute()

println ''
println '|-------------------------------------------------------------------------------------------------------------'
println ''
println "[Info] This is sc_preprocess, latest comitted sha ::: " + latest_sha.text
println '|-------------------------------------------------------------------------------------------------------------'
println ''

ch_fastq    = Channel
                .fromFilePairs(params.fastq, checkIfExists: true)

// Define the final workflow:
workflow SCRNASEQ {

    take:
        genome
        gtf
        fastq

    main:
    
    //-------------------------------------------------------------------------------------------------------------------------------//
    // Create the spliced/unspliced expanded transcriptome:
    include {   ParseExonIntronTx }     from './modules/ParseExonIntronTx'  addParams(  gene_name:  params.ref_gene_name,
                                                                                        gene_id:    params.ref_gene_id,
                                                                                        gene_type:  params.ref_gene_type,
                                                                                        readlength: params.cdna_readlength,
                                                                                        chrM:       params.ref_chrM,
                                                                                        outdir:     params.idx_outdir,
                                                                                        threads:    params.parse_intron_threads,
                                                                                        mem:        params.parse_intron_mem)
    
    // provide reference genome and transcriptome:
    ParseExonIntronTx (genome, gtf) 
                
    //-------------------------------------------------------------------------------------------------------------------------------//
    // Index it:
    include {   AlevinIndex }           from './modules/alevin_index'       addParams(  outdir:     params.idx_outdir,
                                                                                        additional: params.idx_salmonargs,
                                                                                        threads:    params.idx_threads,
                                                                                        mem:        params.idx_mem)

    AlevinIndex(genome,                             // genome fasta file channel
                ParseExonIntronTx.out.txtome,       // expanded transcriptone (exon+intron)
                params.idx_name)                    // output index name

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Quantify with Alevin:
    if(!params.skip_quant){

        include {   AlevinQuant }           from './modules/alevin_quant'       addParams(  outdir:     params.quant_outdir,
                                                                                            libtype:    params.quant_libtype,
                                                                                            additional: params.quant_additional,
                                                                                            threads:    params.quant_threads,
                                                                                            mem:        params.quant_mem)
        
        AlevinQuant(fastq,                              // fastq channel
                    AlevinIndex.out.idx,                // the index itself
                    ParseExonIntronTx.out.tgmap,        // transcript2gene map for gene-level aggregation
                    ParseExonIntronTx.out.rrnagenes,    // list of rRNA gene names
                    ParseExonIntronTx.out.mtgenes)      // list of mito genes

    } else params.skip_mtx = true

    //-------------------------------------------------------------------------------------------------------------------------------//
    if(!params.skip_mtx){
    
        // Load into R, split into spliced/unspliced counts, save as generic MatrixMarket file (gzipped):
        include {   WriteMtx }               from './modules/write_mtx'         addParams(  outdir:     params.mtx_outdir,
                                                                                            threads:    params.mtx_threads,
                                                                                            mem:        params.mtx_mem)
                                                                                            
        WriteMtx(   AlevinQuant.out.quants.map { file -> tuple(file.baseName, file) },
                    ParseExonIntronTx.out.features,     
                    ParseExonIntronTx.out.gene2type)   
    }
   
}                

// Run it:
workflow { 
    SCRNASEQ( params.ref_genome, params.ref_gtf, ch_fastq ) 
}
