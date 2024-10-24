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


process link_transcriptome_novel
{
    publishDir "${params.dump}", mode: 'copy', overwrite: params.force
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
            path("*_complete_txome.json"),
            path("novel_complete_digest/", type: 'dir')
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Writing transcriptome hash to digest"
                echo "Annotation: !{annotation}"
                echo "Transcripts: !{transcripts}"
                echo "Index: !{index}"
            fi
            verbose=""
            if [[ "!{params.log}" == "DEBUG" ]]; then
                verbose="--verbose"
                set -x
            fi

            mkdir -p novel_complete_digest
            
            gzip -cd !{transcripts} \
            > transcripts.fa

            compute_fasta_digest --reference transcripts.fa \
                --out digest.json
            
            sed -e's/seq_hash/SeqHash/' -e's/name_hash/NameHash/' \
                < digest.json \
                > novel_complete_digest/info.json

            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Generating linked transcriptome and txdb for tximeta"
            fi
            Rscript ${verbose} !{projectDir}/bin/R/link_transcriptome.R \
                -a !{annotation} -t !{transcripts} -i novel_complete_digest
        '''
}