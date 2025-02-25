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
 

process correct_flair_annotation
{
    publishDir "${params.dump}/", mode: 'copy', enable: params.dump!='', overwrite: params.force
    if (params.manage_resources)
    {
        cpus 1
        memory '16.GB'
    }
    input:
        path(annotation)
    output:
        path("*_complete_annotation.gtf.gz")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Fixing novel annotation..."
                echo "Annotation: !{annotation}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            python3 !{projectDir}/bin/python/correct_flair_annotation.py !{annotation} \
            | gzip -c \
            > novel_complete_annotation.gtf.gz
        '''
}