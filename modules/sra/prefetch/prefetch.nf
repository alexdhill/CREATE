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


process prefetch
{
    if (params.manage_resources)
    {
        cpus 1
        memory '2.GB'
    }
    errorStrategy "finish"
    maxRetries 10
    input:
        val(acc)
    output:
        tuple(
            val("${acc}"),
            path("${acc}")
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Downloading accession..."
                echo "Sample: !{acc}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            prefetch !{acc} --force ALL --max-size u --output-directory !{acc}
        '''
}
