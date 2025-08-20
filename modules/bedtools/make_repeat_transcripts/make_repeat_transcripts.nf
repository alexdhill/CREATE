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


process make_repeat_transcripts
{
    publishDir "${params.outdir}/", mode: 'copy', enabled: params.keep, overwrite: params.force
    if (params.manage_resources)
    {
        cpus 8
        memory '16.GB' // TODO
    }
    input:
        tuple(
            path(reference),
            path(regions)
        )
    output:
        path("*_repeat_transcripts.fa.gz")
    shell:
        if (params.genome=="T2T")
        {
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Creating repeat transcripts..."
                echo "Repeats: !{regions}"
                echo "Reference: !{reference}"
                echo "Version: 2"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            # Convert BED name
            mkfifo genome regions
            cat !{regions} \
            | awk '{print $1"\t"$2"\t"$3"\t"$4"_range="$1":"$2"-"$3"_strand="$6"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' \
            > regions &
            pigz -cdp !{task.cpus} !{reference} > genome &
            bedtools getfasta -nameOnly -fi genome -bed regions.bed \
            | gzip --best \
            > !{params.genome}v2_repeat_transcripts.fa.gz
        '''
        }
        else
        {
            '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Creating repeat transcripts..."
                echo "Repeats: !{regions}"
                echo "Reference: !{reference}"
                echo "Version: !{params.version}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            mkfifo genome
            pigz -cdp !{task.cpus} !{reference} > genome &
            bedtools getfasta -nameOnly -fi genome -bed !{regions} \
            | pigz -cp !{task.cpus} \
            > !{params.genome}v!{params.version}_repeat_transcripts.fa.gz
        '''
        }
}