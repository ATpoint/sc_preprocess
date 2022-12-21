#! /usr/bin/env nextflow

nextflow.enable.dsl=2

//------------------------------------------------------------------------
// Intro message
//------------------------------------------------------------------------

Date date = new Date()
String datePart = date.format("yyyy-dd-MM -- ")
String timePart = date.format("HH:mm:ss")
def start_date = datePart + timePart

println ""
println "\u001B[33m========================================================================================================================="
println "Pipeline:      sc_preprocess"
println "GitHub:        https://github.com/ATpoint/sc_preprocess"
println "Documentation: https://github.com/ATpoint/sc_preprocess/README.md"
println "Author:        Alexander Toenges (@ATpoint)"
println "Runname:       $workflow.runName"
println "Profile:       $workflow.profile"
println "Start:         $start_date"
println "=========================================================================================================================\u001B[0m"

//------------------------------------------------------------------------
// Validate input params via schema.nf
//------------------------------------------------------------------------

 evaluate(new File("${baseDir}/functions/validate_schema_params.nf"))

//------------------------------------------------------------------------
// Load the modules and pass params
//------------------------------------------------------------------------

include { ParseExonIntronTx } from './modules/parse_exon_intron_tx.nf' addParams(gene_id: params.gene_id,
                                                                                 gene_name: params.gene_name,
                                                                                 gene_type: params.gene_type,
                                                                                 read_length: params.read_length,
                                                                                 chrM: params.chrM,
                                                                                 rrna: params.rrna,
                                                                                 outdir: params.idx_outdir)
                                                                                            
//------------------------------------------------------------------------                                                                                            

include { AlevinIndex } from './modules/alevin_index' addParams(outdir: params.idx_outdir, additional: params.idx_args)

//------------------------------------------------------------------------               

include {ValidateSamplesheet } from './modules/validate_samplesheet.nf'

//------------------------------------------------------------------------          

include { Quant } from './modules/alevin_quant' addParams(outdir: params.quant_outdir,     
                                                          r1_type: params.r1_type,                                   
                                                          r2_type: params.r2_type,
                                                          args: params.quant_args)

//------------------------------------------------------------------------          

include { AlevinIndexFB } from './modules/alevin_index' addParams(outdir: params.idx_outdir)

//------------------------------------------------------------------------

include { QuantFB } from './modules/alevin_quant' addParams(outdir: params.quant_outdir,     
                                                            r1_type: params.r1_type,                                   
                                                            r2_type: params.r2_type_fb,
                                                            args: params.quant_args_fb)

//------------------------------------------------------------------------

include { Mtx } from './modules/write_mtx' addParams(outdir: params.mtx_outdir)

//------------------------------------------------------------------------

include { AlevinQC } from './modules/alevin_qc' addParams(outdir: params.qc_outdir)

//------------------------------------------------------------------------

include { WriteNcells } from './modules/summary_cells' addParams(outdir: params.qc_outdir)

//------------------------------------------------------------------------

include{ CommandLines } from './modules/commandline' addParams(outdir: params.pipe_dir)                                                              

//------------------------------------------------------------------------                                                                                            

workflow INDEXING {

    take:
        genome
        gtf

    main:        
        ParseExonIntronTx (genome, gtf) 
        AlevinIndex(genome, ParseExonIntronTx.out.txtome, params.idx_name)     

    emit:
        idx   = AlevinIndex.out.idx      
        tgmap = ParseExonIntronTx.out.tgmap
        rrnagenes  = ParseExonIntronTx.out.rrnagenes
        mtgenes = ParseExonIntronTx.out.mtgenes
        features  = ParseExonIntronTx.out.features
        gene2type   = ParseExonIntronTx.out.gene2type
        versions = ParseExonIntronTx.out.versions.concat(AlevinIndex.out.versions)

}

workflow VALIDATE {

    take: 
        samplesheet_unvalidated

    main:
        ValidateSamplesheet(samplesheet_unvalidated)

    tuple_fastq = ValidateSamplesheet.out.samplesheet
           .splitCsv(header:true)
           .map {
               
                // Samplesheet allows Nextflow variables to be used, replace by absolute path
                r1 = it['r1']
                        .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                        .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                        .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/"))

                rx = it['r2']
                        .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                        .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                        .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/"))

                r2 = rx.toString()=='' ? '.' : rx   

                // Make sure fastq file paths exist
                is_error = false

                r1_file = file(r1).exists()           
                if(!r1_file){
                    ErrorMessenger("$r1 does not exist")
                    is_error = true
                }

                if(r2!=".") { 
                    r2_file = file(r2).exists()
                    if(!r1_file){
                        ErrorMessenger("$r2 does not exist") 
                        is_error = true
                    }
                }

                // Libtype and whether the fastq file is a gne expression or feature barcode readout
                libtype = it['libtype']
                is_fb = it['is_fb']

                // meta map inspired by nf-core
                meta = [id:it['sample'], libtype: libtype, is_fb: is_fb]     
                reads = [r1: r1, r2: r2]      

                if(!is_error){
                    return tuple(meta, reads)
                } else {
                    return null
                }
                
            }
            .groupTuple(by:0)

    emit:
        tuple_fastq = tuple_fastq
        versions = ValidateSamplesheet.out.versions

}

workflow {

    //--------------------------------------------------------------------
    // Indexing of transcriptome
    //--------------------------------------------------------------------

    if(params.idx==''){
        
        INDEXING(params.genome, params.gtf)
        indexing_versions = INDEXING.out.versions

        indexing_idx       = INDEXING.out.idx      
        indexing_tgmap     = INDEXING.out.tgmap
        indexing_rrnagenes = INDEXING.out.rrnagenes
        indexing_mtgenes   = INDEXING.out.mtgenes
        indexing_features  = INDEXING.out.features
        indexing_gene2type = INDEXING.out.gene2type

    } else {
        
        indexing_versions  = Channel.empty()
        indexing_idx       = params.idx      
        indexing_tgmap     = params.tgmap
        indexing_rrnagenes = params.rrnagenes
        indexing_mtgenes   = params.mtgenes
        indexing_features  = params.features
        indexing_gene2type = params.gene2type

    }

    if(!params.only_idx){

        //--------------------------------------------------------------------
        // Samplesheet validation
        //--------------------------------------------------------------------

        VALIDATE(params.samplesheet)
        validate_versions = VALIDATE.out.versions
    
        ch_geneExpression = VALIDATE.out.tuple_fastq.map{

            aa = it[0]
            r1 = it[1].r1.toList().flatten()
            r2 = it[1].r2.toList().flatten()
            if(aa.is_fb=='false') return [aa, r1, r2] else return null
            
        }

        ch_featureBarcode = VALIDATE.out.tuple_fastq.map{

            aa = it[0]
            r1 = it[1].r1.toList().flatten()
            r2 = it[1].r2.toList().flatten()
            if(aa.is_fb=='true') return [aa, r1, r2] else return null
        }

        //--------------------------------------------------------------------
        // Quantification of gene expression reads
        //--------------------------------------------------------------------

        Quant(ch_geneExpression, indexing_idx.collect(), indexing_tgmap.collect(), indexing_rrnagenes.collect(), indexing_mtgenes.collect())
        quant_versions = Quant.out.versions

        //--------------------------------------------------------------------
        // Indexing of feature barcodes
        //--------------------------------------------------------------------

        if(params.features_file != ''){

            AlevinIndexFB(params.features_file, "idx_fb")
            indexfb_versions = AlevinIndexFB.out.versions

        } else {

            indexfb_versions = Channel.empty()

        }

        //--------------------------------------------------------------------
        // Quantification of feature barcode reads
        //--------------------------------------------------------------------

        if(params.features_file != ''){

            QuantFB(ch_featureBarcode, AlevinIndexFB.out.idx.collect(), AlevinIndexFB.out.tgmap.collect())
            quant_fb_channel = QuantFB.out.quants
            quantfb_versions = QuantFB.out.versions

        } else {

            quant_fb_channel = Channel.empty()
            quantfb_versions = Channel.empty()

        }

        //--------------------------------------------------------------------
        // Splitting of quants into spliced and unspliced counts.
        // If feature barcodes were present for a given sample then 
        // retain only barcodes of the FB that made it to the gene expression
        // matrix. Perform barcode translation if set in the params.
        //--------------------------------------------------------------------

        ch_for_mtx = Quant.out.quants.concat(quant_fb_channel)
                    .map{ [[id: "${it[0].id}"], it[1]] }.groupTuple(by: 0) // retain only .id and the quants channel as part of the tuple

        Mtx(ch_for_mtx, indexing_features.collect(), indexing_gene2type.collect(), params.translate_list)
        mtx_versions = Mtx.out.versions
        
        //--------------------------------------------------------------------
        // AlevinQC reports and cell summary
        //--------------------------------------------------------------------
        
        ch_for_aqc = Quant.out.quants.concat(quant_fb_channel).map{

            ii = it[0].is_fb=='true' ? [[id: "${it[0].id}_fb"], it[1]] : [[id: "${it[0].id}"], it[1]]
            return ii

        }

        AlevinQC(ch_for_aqc)
        alevinqc_versions = AlevinQC.out.versions

        WriteNcells(Mtx.out.ncells.collect())
        ncells_versions = WriteNcells.out.versions

    } else {

        validate_versions = Channel.empty()
        quant_versions    = Channel.empty()
        alevinqc_versions = Channel.empty()
        indexfb_versions  = Channel.empty()
        quantfb_versions  = Channel.empty()
        mtx_versions      = Channel.empty()
        ncells_versions   = Channel.empty()

    }

    //--------------------------------------------------------------------
    // Command lines and software versions
    //--------------------------------------------------------------------

    x_commands = validate_versions.concat(indexing_versions, quant_versions, alevinqc_versions, indexfb_versions, quantfb_versions,
                                          mtx_versions, ncells_versions)
                 .map {it [1]}.flatten().collect()

    x_versions = validate_versions.concat(indexing_versions.first(),
                                          quant_versions.first(), 
                                          alevinqc_versions.first(), 
                                          indexfb_versions, 
                                          quantfb_versions.first(),
                                          mtx_versions.first(), 
                                          ncells_versions)              
                     .map {it [0]}.flatten().collect()

    CommandLines(x_commands, x_versions)

}