process AlevinIndex {

    label 'process_idx'

    errorStrategy 'finish'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

    input:
    path(genome) 
    path(txtome) 
    val(idxname)
        
    output:
    path("decoynames.txt")
    path("gentrome.fa.gz")
    path(idxname), emit: idx
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script: 

    def decoynames  = "decoynames.txt"
    def gentrome    = "gentrome.fa.gz"

    """
    gzip -cd $genome | grep '^>' | cut -d " " -f 1 | awk '{gsub(\">\",\"\");print}' > $decoynames

    cat $txtome $genome > $gentrome

    salmon index --no-version-check \
        -t $gentrome \
        -d $decoynames \
        -i $idxname \
        -p $task.cpus \
        $params.additional

    echo ${task.process}: > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'salmon:' \$(salmon --version | cut -d " " -f2) > versions.txt        
    """                

}

process AlevinIndexFB {

    cpus 1
    memory 500.MB

    errorStrategy 'finish'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

    input:
    path(featureset) 
    val(idxname)
        
    output:
    path(idxname), emit: idx
    path("$idxname/tgmap.txt"), emit: tgmap
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script: 
    """
    salmon index --no-version-check -t ${featureset} -i ${idxname} --features -k7

    awk 'OFS=\"\t\" {print \$1, \$1}' $featureset > "$idxname/tgmap.txt"

    echo ${task.process}: > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'salmon:' \$(salmon --version | cut -d " " -f2) > versions.txt
    """                

}