/*
 * REQUIRED NOTICE: Copyright (c) 2020-2023, Regents of the University of California
 * All rights reserved. https://polyformproject.org/licenses/noncommercial/1.0.0
 * 
 * This software was developed by the Daniel Kim lab at the University of California, Santa Cruz.
 * Authors: Roman E. Reggiardo, Vikas Peddu, Alex D. Hill
 * 
 * The licensor grants you a copyright license for the software to do everything you might do with
 * the software that would otherwise infringe the licensor’s copyright in it for any permitted
 * purpose.
 * 
 * As far as the law allows, the software comes as is, without any warranty or condition, and the
 * licensor will not be liable to you for any damages arising out of these terms or the use or
 * nature of the software, under any kind of legal claim.
 */
 

process salmon_quant_novel
{
    publishDir "${params.outdir}/quant/", mode: 'copy', enable: params.keep, overwrite: params.force
    errorStrategy 'ignore'
    if (params.manage_resources)
    {
        cpus 8
        memory '32.GB'
        time '6h'
    }
    input:
        tuple(
            val(sample),
            path(read_1),
            path(read_2),
            val(nreads),
            path(index)
        )
    output:
        path("${sample}/")
    shell:
        '''
            if [[ "!{params.log}"=="INFO" || "!{params.log}"=="DEBUG" ]]; then
                echo "Quantifying paired reads"
                echo "Sample: !{sample} (!{nreads} reads)"
                echo "Read 1: !{read_1}"
                echo "Read 2: !{read_2}"
                echo "Index: !{index}"
            fi
            if [[ !{params.log}=="DEBUG" ]]; then
                set -x
            fi

            salmon quant --libType A -1 !{read_1} -2 !{read_2} \
                -i !{index} -p 8 --output !{sample} \
                --seqBias --gcBias --writeUnmappedNames \
                --validateMappings --recoverOrphans --rangeFactorizationBins 4
            
            if [[ ! -e !{sample}/quant.sf ]]; then
                echo "\033[1;31mERR: Salmon quantification failed\033[0m" 1>&2
                exit 1
            fi
        '''
}