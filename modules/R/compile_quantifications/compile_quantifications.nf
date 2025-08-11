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


process compile_quantifications
{
    publishDir "${params.outdir}/", mode: 'copy', overwrite: params.force
    if (params.manage_resources)
    {
        cpus 1
        memory '24GB'
    }
    input:
        tuple(
            path(quants),
            path(reference),
            path(metadata)
        )
    output:
        path("*counts/")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Compiling quants to H5 SummarizedExperiment"
                echo "Reference: !{reference}"
                echo "Metadata: !{metadata}"
                echo "Quants:"
                sed 's/ /\\n/g' <<< "!{quants}"
            fi
            verbose=""
            if [[ "!{params.log}" == "DEBUG" ]]; then
                verbose="--verbose"
                set -x
            fi

            splintr=""
            if [ -n "!{params.get('library')}" ] && [ "!{params.library}" == "single_cell" ]; then
                splintr="-s"
            fi

            transcripts=""
            if [ -n "!{params.get('export_transcripts')}" ]; then
                transcripts="-t"
            fi

            mkdir -p quants && mv !{quants} quants/
            Rscript ${verbose} !{projectDir}/bin/R/compile_quantifications.R \
                -q quants -r !{reference} -m !{metadata} ${splintr} ${transcripts}
        '''
}