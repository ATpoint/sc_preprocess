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

include { ValidateSamplesheet   }   from './modules/validate_samplesheet.nf'    addParams(  outdir: params.outdir)

//------------------------------------------------------------------------

include { ParseExonIntronTx     }   from './modules/parse_exon_intron_tx.nf'    addParams(  gene_name:      params.gene_name,
                                                                                            gene_id:        params.gene_id,
                                                                                            gene_type:      params.gene_type,
                                                                                            chrM:           params.chrM,
                                                                                            rrna:           params.rrna,
                                                                                            outdir:         params.idx_outdir)
                                                                                            
//------------------------------------------------------------------------                                                                                            

include { AlevinIndex           }   from './modules/alevin_index'               addParams(  outdir:         params.idx_outdir,
                                                                                            additional:     params.idx_args)

//------------------------------------------------------------------------

include { AlevinIndexFB         }   from './modules/alevin_index'               addParams(  outdir:         params.idx_outdir)

//------------------------------------------------------------------------

// combine the read geometry and additional params into a single string
use_quant_args = params.quant_args + ' ' + params.r1_type + ' ' + params.r2_type
if(use_quant_args == ' ') use_quant_args = ''

include { AlevinQuant           }   from './modules/alevin_quant'               addParams(  outdir:         params.quant_outdir,
                                                                                            libtype:        params.libtype,
                                                                                            additional:     use_quant_args)

//------------------------------------------------------------------------

// combine the read geometry and additional params into a single string
use_quant_args_sf = params.quant_args + ' ' + params.r1_type + ' ' + params.r2_type_fb
if(use_quant_args_sf == ' ') use_quant_args_sf = ''

include { AlevinQuantFB         }   from './modules/alevin_quant'               addParams(  outdir:         params.quant_outdir,
                                                                                            libtype:        params.libtype,  
                                                                                            suffix:         params.fb_suffix,
                                                                                            additional:     use_quant_args_sf)

//------------------------------------------------------------------------                                                                                            

include { WriteMtx              }   from './modules/write_mtx'                  addParams(  outdir:         params.mtx_outdir)

//------------------------------------------------------------------------  

include { WriteNcells           }   from './modules/summary_cells'              addParams(  outdir:         params.qc_outdir)

//------------------------------------------------------------------------  

include { AlevinQC              }   from './modules/alevin_qc'                  addParams(  outdir:         params.qc_outdir)
                                                                                            
//------------------------------------------------------------------------      
// Validate that samplesheet exists and that the fastq files exist
//------------------------------------------------------------------------

if(!params.idx_only){

    // Validate that fastq files in samplesheet exist as files on disk

    fastq_no_exist  = [:]
    
    if(!(new File(params.samplesheet)).exists()){
        println "\u001B[31m========================================================================================================================="
        println "[VALIDATION ERROR]"
        println "The samplesheet does not exist!"
        println "=========================================================================================================================\u001B[0m"
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
        println "\u001B[31m========================================================================================================================="
        println "[VALIDATION ERROR]"
        println "The following fastq files in the samplesheet do not exist:"
        println fastq_no_exist
        println "=========================================================================================================================\u001B[0m"
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
            println "\u001B[31m========================================================================================================================="
            println "[VALIDATION ERROR]"
            println "${nmm} does not exist!"
            println "=========================================================================================================================\u001B[0m"
            
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
workflow INDEX_FB {

    take:
        fb_file
        idxname

    main:
        AlevinIndexFB(fb_file, idxname)

    emit: 
        idx    = AlevinIndexFB.out.idx
        tgmap  = AlevinIndexFB.out.tgmap

}

// Quant against gentrome
workflow QUANT {

    take: 
        fastq
        tgmap  
        rrnagenes
        mtgenes
        idx

    main:
        AlevinQuant(fastq, idx, tgmap, rrnagenes, mtgenes)

    emit:
        quants = AlevinQuant.out.quants     

}

// Quant against feature barcodes
workflow QUANT_FB {

    take: 
        fastq  
        idx
        tgmap

    main:
        AlevinQuantFB(fastq, idx, tgmap)

    emit:
        quants = AlevinQuantFB.out.quants                

}

// Split quants into spliced/unspliced and write as mtx
workflow WRITE_MTX {

    take:
        quants
        features
        gene2type
        samplesheet

    main:    
        WriteMtx(quants, features, gene2type, samplesheet) 

    emit:
        mtx      = WriteMtx.out.mtx
        barcodes = WriteMtx.out.barcodes      
        features = WriteMtx.out.features
        ncells   = WriteMtx.out.ncells

}

workflow SUMMARY {
    
    take:
        dirs

    main:
        WriteNcells(dirs)
            
}

// Basic QC using alevinQC (Bioc)
workflow ALEVIN_QC {

    take:
        quants

        main:
            AlevinQC(quants, 'empty')

}

workflow ALEVIN_QC_FB {

    take:
        quants

        main:
            AlevinQC(quants, params.fb_suffix)

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
            if(k['is_fb']=='false') {
                tuple(k['sample_id'], k['R1'], k['R2'], k['is_fb'])
            } else null
        }.groupTuple(by: 0)

        // channel with the feature barcode fastq pairs:
        ch_input_quant_fb = ch_samplesheet.map { k -> 
            if(k['is_fb']=="true") {
                tuple(k['sample_id'], k['R1'], k['R2'], k['is_fb'])
            } else null
        }.groupTuple(by: 0)

        // Index gentrome and then quantify against it:
        QUANT(ch_input_quant, tgmap, rrna, mtrna, idx)

        ALEVIN_QC(QUANT.out.quants)

        // Index SFs and then quantify against it:
        if(params.features_file!=''){

            INDEX_FB(params.features_file, "idx_features")

            QUANT_FB(ch_input_quant_fb, INDEX_FB.out.idx, INDEX_FB.out.tgmap)

            ALEVIN_QC_FB(QUANT_FB.out.quants)
            
        }

        /* 
         * Split alevin quantifications into spliced/unspliced.
         * optionally match these with the feature barcode libraries and retain only those CBs present in both RNA and FBs,
         * optionally perform barcode translation,
         * and save as mtx.gz and the barcodes/features as tsv.gz
         * For this combine the RNA and FB alevin quants into a single tuple
         */

        // translation table, only used if translate_barcodes is true
        if(params.translate_barcodes==true){
                ch_translate = Channel.fromPath(params.translate_list, checkIfExists: true)
        } else ch_translate = Channel.fromPath("/")

        // combine RNA and FB alevins into a tuple:
        if(params.features_file==''){
            ch_quants = QUANT.out.quants.map { qm -> tuple(qm[0], qm[1], "/")}
        } else {
            ch_quants = QUANT.out.quants.join(QUANT_FB.out.quants, remainder: true).map { boom ->
                if(boom[2]==null) {
                    return tuple(boom[0], boom[1], "/")
                } else return boom
            }
        }

        WRITE_MTX(ch_quants, ftrs, gene2type, ch_translate.collect())
        WRITE_MTX.out.ncells.collect().view()
        SUMMARY(WRITE_MTX.out.ncells.collect())

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
        println "\u001B[33m========================================================================================================================="
        println "Pipeline completed!"
        println "End: $end_date"
        println "Results are in:"
        println od
        //println "A summary file with cellnumbers per sample is at:"
        //println "$smyfile"
        println "=========================================================================================================================\u001B[0m"
        println ""
    }

}