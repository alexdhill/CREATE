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
 

process make_transcript_map
{
    publishDir "${params.outdir}/", mode: 'copy', overwrite: params.force
    input:
        tuple(
            path(gencode_annotation),
            path(repeat_annotation)
        )
    output:
        path("*_map.tx2g")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Generating transcript to gene map..."
                echo "Repeat Annotation: !{repeat_annotation}"
                echo "Gencode Annotation: !{gencode_annotation}"
            fi
            version="!{params.version}"
            if [[ "!{params.genome}"=="T2T" ]]; then
                version="2"
            fi
            verbose=""
            if [[ "!{params.log}" == "DEBUG" ]]; then
                verbose="--verbose"
                set -x
            fi

            gzip -cd !{gencode_annotation} !{repeat_annotation} \
            | grep -v "^#" \
            | awk -F'\\t' '{print $3"\\t"$9}' \
            | grep '^transcript' \
            | awk -F'\\t' '{print $2}' \
            | sed 's/gene_type/gene_biotype/' \
            > transcript_info.txt
            Rscript ${verbose} !{projectDir}/bin/R/parse_transcript_info.R transcript_info.txt !{params.genome}v${version}_complete_map.tx2g
        '''
}