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
 

process salmon_index
{
    publishDir "${params.outdir}/", mode: 'copy', overwrite: params.force
    if (params.manage_resources)
    {
        cpus 8
        memory '24.GB'
    }
    input:
        tuple(
            path(transcripts),
            val(prefix),
        )
    output:
        path("*.sidx")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Creating complete salmon index"
                echo "Transcripts: !{transcripts}"
            fi
            version="!{params.version}"
            if [[ "!{params.genome}"=="T2T" ]]; then
                version="2"
            fi
            if [[ "!{params.log}"=="DEBUG" ]]; then
                set -x
            fi

            zcat !{transcripts} > txome.fa
            salmon index -t txome.fa \
                -p 8 \
                -i !{params.genome}v${version}_!{prefix}_index_v$(salmon --version | awk '{print $2}').sidx
        '''
}