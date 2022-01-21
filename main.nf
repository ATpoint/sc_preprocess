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
println "\u001B[33m======================================================================"
println "Pipeline:  sc_preprocess"
println "GitHub:    https://github.com/ATpoint/sc_preprocess"
println "Author:    Alexander Toenges (@ATpoint)"
println "Runname:   $workflow.runName"
println "Profile:   $workflow.profile"
println "Start:     $start_date"
println "======================================================================\u001B[0m"

//------------------------------------------------------------------------
// Validate input params via schema.nf
//------------------------------------------------------------------------

 evaluate(new File("${baseDir}/functions/validate_schema_params.nf"))

//------------------------------------------------------------------------
// Load the modules and pass params
//------------------------------------------------------------------------

include { ValidateSamplesheet   }   from './modules/validate_samplesheet.nf'    addParams(  outdir: params.outdir)

include { ParseExonIntronTx     }   from './modules/parse_exon_intron_tx.nf'    addParams(  gene_name:      params.gene_name,
                                                                                            gene_id:        params.gene_id,
                                                                                            gene_type:      params.gene_type,
                                                                                            chrM:           params.chrM,
                                                                                            rrna:           params.rrna,
                                                                                            outdir:         params.idx_outdir)

include { AlevinIndex           }   from './modules/alevin_index'               addParams(  outdir:         params.idx_outdir,
                                                                                            additional:     params.idx_args)

include { AlevinIndexSF         }   from './modules/alevin_index'               addParams(  outdir:         params.idx_outdir)



include { AlevinQuant           }   from './modules/alevin_quant'               addParams(  outdir:         params.quant_outdir,
                                                                                            libtype:        params.quant_libtype,
                                                                                            additional:     params.quant_args)

include { AlevinQuantSF         }   from './modules/alevin_quant'               addParams(  outdir:         params.quant_outdir,
                                                                                            libtype:        params.quant_libtype,  
                                                                                            suffix:         params.quant_sf_suffix,
                                                                                            additional:     params.quant_sf_args)

include { WriteMtx              }   from './modules/write_mtx'                  addParams(  outdir:         params.mtx_outdir)

include { WriteMtxSF            }   from './modules/write_mtx'                  addParams(  outdir:         params.mtx_outdir,
                                                                                            suffix:         params.quant_sf_suffix)

include { AlevinQC              }   from './modules/alevin_qc'                  addParams(  outdir:         params.qc_outdir)
                                                                                            
//------------------------------------------------------------------------      
// Validate that samplesheet exists and that the fastq files exist
//------------------------------------------------------------------------

if(!params.idx_only){

    // Validate that fastq files in samplesheet exist as files on disk

    fastq_no_exist  = [:]
    
    if(!(new File(params.samplesheet)).exists()){
        println "\u001B[31m======================================================================"
        println "[VALIDATION ERROR]"
        println "The samplesheet does not exist!"
        println "======================================================================\u001B[0m"
        System.exit(1)
    }

    def file_samplesheet = new FileReader(params.samplesheet)

    file_samplesheet.eachLine(1) {line, number-> 

        if(number>1){

            s = line
                    .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                    .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                    .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/")).split(',')
                
            if(!new File(s[1]).exists()) fastq_no_exist.put(s[0], s[1])
            if(!new File(s[2]).exists()) fastq_no_exist.put(s[0], s[2])
                    
        }

    }

    if(fastq_no_exist.size() > 0){
        println "\u001B[31m======================================================================"
        println "[VALIDATION ERROR]"
        println "The following fastq files in the samplesheet do not exist:"
        println fastq_no_exist
        println "======================================================================\u001B[0m"
        System.exit(1)
    }

    // Channel to samplesheet:
    if(params.samplesheet != '') { 
        ch_samplesheet_in = Channel.fromPath(params.samplesheet, checkIfExists: true) 
    }

}

//------------------------------------------------------------------------      
// Validate that (in case a premade index is provided) the individual
// elements exist
//------------------------------------------------------------------------

 def ConvertBool2String(indata='') {
    if(indata instanceof Boolean){
        return ''
    } else {
        return indata
    }
 }

if(params.idx != ''){

        // little hack when the param is empty (which groovy then converts to boolean), convert to string
        // so File() will not complain about Boolean rather than having String input
        use_idx   = ConvertBool2String(params.idx)
        use_tgmap = ConvertBool2String(params.tgmap)
        use_rrna  = ConvertBool2String(params.rrnagenes)
        use_mtrna = ConvertBool2String(params.mtrnagenes)
        use_expanded_features = ConvertBool2String(params.expanded_features)
        use_gene2type = ConvertBool2String(params.gene2type)

        // validate existance:
        def not_exist = [:]
        if(!new File(use_idx).exists()) not_exist.put("use_idx", use_idx)
        not_exist.each { mm, nn -> 
        
            def nmm = mm.replaceAll("use_", "--")
            println "\u001B[31m======================================================================"
            println "[VALIDATION ERROR]"
            println "${nmm} does not exist!"
            println "======================================================================\u001B[0m"
            
        }

        if(not_exist.size() > 0) System.exit(1)

    }

//------------------------------------------------------------------------      
// Define subworkflows
//------------------------------------------------------------------------

// Samplesheet validation
workflow VALIDATE {

    take:
        samplesheet
        
    main:        
        ValidateSamplesheet(samplesheet, baseDir, launchDir, projectDir)

    emit:
        samplesheet = ValidateSamplesheet.out.ssheet

}

// Indexing of expanded transcriptome + genome (=gentrome)
workflow INDEX_GENTROME {

    take:
        genome
        gtf

    main:        
        ParseExonIntronTx (genome, gtf) 
        AlevinIndex(genome, ParseExonIntronTx.out.txtome, params.idx_name)     

     emit:
        tgmap = ParseExonIntronTx.out.tgmap
        rrna  = ParseExonIntronTx.out.rrnagenes
        mtrna = ParseExonIntronTx.out.mtgenes
        idx   = AlevinIndex.out.idx      
        ftrs  = ParseExonIntronTx.out.features
        g2t   = ParseExonIntronTx.out.gene2type

}

// Indexing of feature barcodes
workflow INDEX_SF {

    take:
        sf_file
        idxname

    main:
        AlevinIndexSF(sf_file, idxname)

    emit: 
        idx    = AlevinIndexSF.out.idx
        tgmap  = AlevinIndexSF.out.tgmap

}

// Quant against gentrome
workflow QUANT {

    take: 
        fastq
        tgmap  
        rrnagenes
        mtgenes
        idx_gentrome

    main:
        AlevinQuant(fastq, idx_gentrome, tgmap, rrnagenes, mtgenes)

    emit:
        quants = AlevinQuant.out.quants     

}

// Quant against feature barcodes
workflow QUANT_SF {

    take: 
        fastq  
        idx
        tgmap
        whitelist

    main:
        AlevinQuantSF(fastq, idx, tgmap, whitelist)

    emit:
        quants = AlevinQuantSF.out.quants                

}

// Split quants into spliced/unspliced and write as mtx
workflow WRITE_MTX {

    take:
        quants
        features
        gene2type

    main:    
        WriteMtx(quants, features, gene2type) 

    emit:
        barcodes = WriteMtx.out.barcodes      

}

// Write feature barcode counts as mtx
workflow WRITE_MTX_SF {

    take:
        quants
        suffix

    main:    
        WriteMtxSF(quants, suffix) 

}

// Basic QC using alevinQC (Bioc)
workflow ALEVIN_QC {

    take:
        quants

        main:
            AlevinQC(quants)

}

//------------------------------------------------------------------------
// Assemble main workflow from subworkflows
//------------------------------------------------------------------------

workflow SC_PREPROCESS { 
    
    take:
        samplesheet
        idx
        tgmap
        rrna
        mtrna
        gene2type
        ftrs

    main:

        VALIDATE(samplesheet)

        // channel with the validated samplesheet
        ch_samplesheet = VALIDATE.out.samplesheet.splitCsv(header: true)

        // channel with the transcriptomic fastq pairs:
        ch_input_quant = ch_samplesheet.map { k -> 
            if(k['is_sf']=='false') {
                tuple(k['sample_id'], k['R1'], k['R2'], k['is_sf'])
            } else null
        }.groupTuple(by: 0)

        // channel with the feature barcode fastq pairs:
        ch_input_quant_sf = ch_samplesheet.map { k -> 
            if(k['is_sf']=="true") {
                tuple(k['sample_id'], k['R1'], k['R2'], k['is_sf'])
            } else null
        }.groupTuple(by: 0)

        // Index gentrome and then quantify against it:
        QUANT(ch_input_quant, tgmap, rrna, mtrna, idx)

        WRITE_MTX(QUANT.out.quants, ftrs, gene2type)

        ALEVIN_QC(QUANT.out.quants)

        // Index SFs and then quantify against it:
        if(params.features_file!=''){

            INDEX_SF(params.features_file, "idx_features")

            QUANT_SF(ch_input_quant_sf, INDEX_SF.out.idx, INDEX_SF.out.tgmap, WRITE_MTX.out.barcodes)

            WRITE_MTX_SF(QUANT_SF.out.quants, params.quant_sf_suffix)
            
        }

}

//------------------------------------------------------------------------
// Run the main workflow
//------------------------------------------------------------------------

workflow { 
    
    // Run indexing of gentrome if no premade index is provided:
    if(params.idx == ''){
    
        INDEX_GENTROME(params.genome, params.gtf)
        use_idx       = INDEX_GENTROME.out.idx
        use_tgmap     = INDEX_GENTROME.out.tgmap
        use_rrna      = INDEX_GENTROME.out.rrna
        use_mtrna     = INDEX_GENTROME.out.mtrna
        use_expanded_features = INDEX_GENTROME.out.ftrs   // expanded features
        use_gene2type = INDEX_GENTROME.out.g2t

    }
    
    // Run the quantification, mtx and QC:
    if(!params.idx_only){

        SC_PREPROCESS(params.samplesheet, use_idx, use_tgmap, use_rrna, use_mtrna, use_gene2type, use_expanded_features) 

    }

    // Final message printing the output directory:
    def od = params.outdir
    workflow.onComplete {
        Date date2 = new Date()
        String datePart2 = date2.format("yyyy-dd-MM -- ")
        String timePart2 = date2.format("HH:mm:ss")
        def end_date = datePart2 + timePart2
        println ""
        println "\u001B[33m======================================================================"
        println "Pipeline completed!"
        println "End: $end_date"
        println "Results are in:"
        println od
        println "======================================================================\u001B[0m"
        println ""
    }

}