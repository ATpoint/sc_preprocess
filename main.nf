#! /usr/bin/env nextflow

nextflow.enable.dsl=2

//------------------------------------------------------------------------

evaluate(new File("${baseDir}/functions/validate_schema_params.nf"))

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
                                                                                            suffix:         params.quant_sf_sfx,
                                                                                            additional:     params.quant_sf_args)

include { WriteMtx              }   from './modules/write_mtx'                  addParams(  outdir:         params.mtx_outdir)

include { WriteMtxSF            }   from './modules/write_mtx'                  addParams(  outdir:         params.mtx_outdir)

include { AlevinQC              }   from './modules/alevin_qc'                  addParams(  outdir:         params.qc_outdir)
                                                                                            
//------------------------------------------------------------------------      

// Validate that fastq files in samplesheet exist as files on disk

fastq_no_exist  = [:]

new FileReader(params.samplesheet).eachLine(1) {line, number-> 

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


workflow VALIDATE {

    take:
        samplesheet
        
    main:        
        ValidateSamplesheet(samplesheet, baseDir, launchDir, projectDir)

    emit:
        samplesheet = ValidateSamplesheet.out.ssheet

}

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

workflow QUANT_SF {

    take: 
        fastq  
        idx
        tgmap

    main:
        AlevinQuantSF(fastq, idx, tgmap)

    emit:
        quants = AlevinQuantSF.out.quants                

}

workflow WRITE_MTX {

    take:
        quants
        features
        gene2type

    main:    
        WriteMtx(quants, features, gene2type) 

}

workflow WRITE_MTX_SF {

    take:
        quants

    main:    
        WriteMtxSF(quants) 

}

workflow ALEVIN_QC {

    take:
        quants

        main:
            AlevinQC(quants)

}

//------------------------------------------------------------------------

// samplesheet channel:
if(params.samplesheet != '') { ch_samplesheet_in = Channel.fromPath(params.samplesheet, checkIfExists: true) }

//------------------------------------------------------------------------

// Workflow in a logical order:
workflow SC_PREPROCESS { 
    
    // Before doing anything else, validate the samplesheet:
    if(params.samplesheet != ''){

        VALIDATE(ch_samplesheet_in)

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

    }

    // Index gentrome and then quantify against it:
    if(params.samplesheet != ''){   

        INDEX_GENTROME(params.genome, params.gtf)

        QUANT(ch_input_quant, INDEX_GENTROME.out.tgmap, INDEX_GENTROME.out.rrna, 
              INDEX_GENTROME.out.mtrna, INDEX_GENTROME.out.idx)

        WRITE_MTX(QUANT.out.quants, INDEX_GENTROME.out.ftrs, INDEX_GENTROME.out.g2t)

        //ch_for_alevinqc = QUANT.out.quants.flatMap { y ->
        //    y[1]
        //}.collect()

        ALEVIN_QC(QUANT.out.quants)

    }

    // Index SFs and then quantify against it:
    if(params.features_file!=''){ // & !ch_input_quant_sf.toList().isEmpty()){

        INDEX_SF(params.features_file, "idx_sf")

        QUANT_SF(ch_input_quant_sf, INDEX_SF.out.idx, INDEX_SF.out.tgmap)

        WRITE_MTX_SF(QUANT_SF.out.quants)
        
    }

}

workflow { 
    
    SC_PREPROCESS() 
    
    def od = params.outdir
    workflow.onComplete {
        println "\u001B[32m"
        println "Pipeline completed!"
        println ""
        println "Results are in:"
        println od
        println "\u001B[0m"
        println ""
    }

}