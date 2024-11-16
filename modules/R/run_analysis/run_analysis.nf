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
 

process run_analysis
{
    publishDir "${params.outdir}/results", mode: 'copy', overwrite: params.force
    input:
        tuple(
            path(counts),
            path(metadata)
        )
    output:
        tuple(
            path("plots/", optional: true),
            path("data/", optional: true)
        )
    shell:
        if (params.metadata!="")
        {
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Generating biotypes plots..."
                echo "H5 Directory: !{counts}"
            fi
            verbose=""
            if [[ "!{params.log}" == "DEBUG" ]]; then
                verbose="--verbose"
                set -x
            fi

            Rscript ${verbose} !{projectDir}/bin/R/run_analysis.R \
                --helpers !{projectDir}/bin/R/create_helpers.R \
                --quants !{counts} \
                --output .
        '''
        }
}