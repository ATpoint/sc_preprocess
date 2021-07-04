// Alevin/Salmon indexing module:
process AlevinIndex {

    label params.label

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path(genome) // genome fasta channel
    path(txtome) // expanded (exon+intron) transcriptome from ParseExonIntronTx
    val(idxname)
        
    output:
    path("decoynames.txt")
    path("gentrome.fa.gz")
    path(idxname), emit: idx
    
    script: 

    // nf-core-inspired (salmon module)
    def decoynames  = "decoynames.txt"
    def gentrome    = "gentrome.fa.gz"

    """
    zgrep '^>' $genome | cut -d " " -f 1 | awk '{gsub(\">\",\"\");print}' > $decoynames

    cat $txtome $genome > $gentrome

    salmon index --no-version-check \
        -t $gentrome \
        -d $decoynames \
        -i $idxname \
        -p $task.cpus \
        $params.additional
    """                

}