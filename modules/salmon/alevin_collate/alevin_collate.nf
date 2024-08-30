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
 

process alevin_collate
{
    publishDir "${params.outdir}/", mode: 'copy', enable: params.keep, overwrite: params.force
    if (params.manage_resources)
    {
        cpus 8
    }
    input:
        tuple(
            val(sample),
            path(alevin),
            path(whitelist)
        )
    output:
        tuple(
            val("${sample}"),
            path("fry/${sample}")
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Finding barcodes..."
                echo "Sample: !{sample}"
                echo "Alevin: !{alevin}"
                echo "Barcode whitelist: !{whitelist}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            mkdir -p fry
            alevin-fry generate-permit-list \
                --input !{alevin} \
                --unfiltered-pl !{whitelist} \
                --expected-ori fw \
                --output-dir fry/!{sample}

            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Collating RAD file..."
            fi
            alevin-fry collate \
                -i fry/!{sample} \
                -r !{alevin} \
                -t !{task.cpus}
        '''
}