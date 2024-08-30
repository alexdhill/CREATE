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


process link_transcriptome
{
    publishDir "${params.outdir}", mode: 'copy', overwrite: params.force
    if (params.manage_resources)
    {
        cpus 1
        memory '16.GB'
    }
    input:
        tuple(
            path(annotation),
            path(transcripts),
            path(index)
        )
    output:
        tuple(
            path(".bfccache/"),
            path("${params.genome}v*_complete_digest/"),
            path("${params.genome}v*_complete_txome.json"),
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Writing transcriptome hash to digest"
                echo "Annotation: !{annotation}"
                echo "Transcripts: !{transcripts}"
                echo "Index: !{index}"
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

            gzip -cd !{transcripts} > txome.fa
            mkdir -p !{params.genome}v${version}_complete_digest
            compute_fasta_digest --reference txome.fa \
                --out digest.json
            sed -e's/seq_hash/SeqHash/' -e's/name_hash/NameHash/' \
                < digest.json \
                > !{params.genome}v${version}_complete_digest/info.json

            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Generating linked transcriptome and txdb for tximeta"
            fi
            Rscript ${verbose} !{projectDir}/bin/R/link_transcriptome.R \
                -a !{annotation} -t !{transcripts} -i !{params.genome}v${version}_complete_digest -p !{params.index}
        '''
}