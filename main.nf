#! /usr/bin/env nextflow

/*  
    ---------------------------------------------------------------------------------------------------------------
    Workflow for 10X scRNA-seq:
    --||    Creation of an index with salmon comprising the spliced- and unspliced transcriptome
            with the entire genome as decoy, starting from GENCODE reference files.

    --||    Quantification of the fastq files with Salmon/Alevin.

    --||    Read into R using tximeta, create spliced- and unspliced count tables, save as mtx.gz for downstream
            analysis.

    --||    Perform basic QC on the data and create a Rmarkdown summary report.                
    ---------------------------------------------------------------------------------------------------------------
*/ 

nextflow.enable.dsl=2

ch_fastq    = Channel
                .fromFilePairs(params.fastq, checkIfExists: true)
                .ifEmpty("No fastq files found")   

// Define the final workflow:
workflow SCRNASEQ {

    take:
        genome
        txtome
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
                                                                                        outdir:     params.idx_outdir)
    
    // provide reference genome and transcriptome:
    ParseExonIntronTx (genome, gtf) 
                
    //-------------------------------------------------------------------------------------------------------------------------------//
    // Index it:
    include {   AlevinIndex }           from './modules/alevin_index'       addParams(  threads:    params.idx_threads,
                                                                                        outdir:     params.idx_outdir,
                                                                                        additional: params.idx_salmonargs,
                                                                                        label:      params.idx_label)

    AlevinIndex(genome,                          // genome fasta file channel
                ParseExonIntronTx.out.txtome,       // expanded transcriptone (exon+intron)
                params.idx_name)                    // output index name

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Quantify with Alevin:
    include {   AlevinQuant }           from './modules/alevin_quant'       addParams(  outdir:     params.quant_outdir,
                                                                                        threads:    params.quant_threads,
                                                                                        libtype:    params.quant_libtype,
                                                                                        additional: params.quant_additional,
                                                                                        label:      params.quant_label)
    
    AlevinQuant(fastq,                                        // fastq channel
                AlevinIndex.out.idx,                // the index itself
                ParseExonIntronTx.out.tgmap,        // transcript2gene map for gene-level aggregation
                ParseExonIntronTx.out.rrnagenes,    // list of rRNA gene names
                ParseExonIntronTx.out.mtgenes)      // list of mito genes

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Load into R, split into spliced/unspliced, save as mtx.gz:
    include {   WriteMtx }               from './modules/write_mtx'         addParams(  outdir:     params.mtx_outdir,
                                                                                        author:     params.author)
                                                                                        
    WriteMtx(   AlevinQuant.out.quants.map { file -> tuple(file.baseName, file) },
                ParseExonIntronTx.out.features,     
                ParseExonIntronTx.out.gene2type)   
   
}                

// Run it:
workflow { 
    SCRNASEQ( params.ref_genome, params.ref_txtome, params.ref_gtf, ch_fastq ) 
}
