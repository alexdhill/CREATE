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
 

process make_novel_reference
{
    publishDir "${params.dump}/", mode: 'copy', enable: params.dump!='', overwrite: params.force
    if (params.manage_resources)
    {
        cpus 1
        memory '1.GB'
    }
    input:
        tuple(
            path(transcripts),
            path(regions),
            path(annotation),
            path(read_map),
            path(reference)
        )
    output:
        tuple(
            path("${params.genome}v${params.genome=='T2T'?'2':params.version}_complete_regions.bed.gz"),
            path("${params.genome}v${params.genome=='T2T'?'2':params.version}_complete_readmap.txt"),
            path("${params.genome}v${params.genome=='T2T'?'2':params.version}_genome.fa.gz")
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Generating novel annotation..."
                echo "Transcripts: !{transcripts}"
                echo "Regions: !{regions}"
                echo "Annotation: !{annotation}"
                echo "Read map: !{read_map}"
                echo "Reference: !{reference}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            mkdir novel_!{reference}
            gzip -c !{regions} > novel_complete_regions.bed.gz
            cp !{read_map} novel_complete_readmap.txt
            cp !{reference}/*genome.fa.gz .
        '''
}