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


process salmon_quant_nanopore
{
    publishDir "${params.outdir}/quant/", mode: 'copy', enable: params.keep, overwrite: params.force
    if (params.manage_resources)
    {
        cpus 8
        memory '32.GB'
        time '6h'
    }
    input:
        tuple(
            val(sample),
            val(nreads),
            path(bam),
            path(reference)
        )
    output:
        path("${sample}/")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Quantifying nanopore reads"
                echo "Sample: !{sample} (!{nreads} reads)"
                echo "BAM: !{bam}"
                echo "Reference: !{reference}"
                echo "User parameters: $(jq '.salmon' !{params.parameters})"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi
            params="$(jq '.salmon' !{params.parameters})"
            if [[ "${params}" == "null" ]]; then
                params=""
            fi

            salmon quant --libType U --ont -a !{bam} \
                -p 8 --noLengthCorrection --noErrorModel \
                -t !{reference}/*complete_transcripts.fa.gz --output !{sample} \
                ${params}

            if [[ ! -e !{sample}/quant.sf ]]; then
                echo "\033[1;31mERR: Salmon quantification failed\033[0m" 1>&2
                exit 1
            fi

            grep 'SeqHash' !{reference}/*complete_digest/info.json \
            | awk '{print $2}' \
            | sed -e's/\"//g' \
            | xargs -I'%' \
                sed -i \
                    's/\"index_seq_hash\": \"\"/\"index_seq_hash\": \"%\"/' \
                    !{sample}/aux_info/meta_info.json

            grep 'NameHash' !{reference}/*complete_digest/info.json \
            | awk '{print $2}' \
            | sed -e's/\"//g' \
            | xargs -I'%' \
                sed -i \
                    's/\"index_name_hash\": \"\"/\"index_name_hash\": \"%\"/' \
                    !{sample}/aux_info/meta_info.json
        '''
}