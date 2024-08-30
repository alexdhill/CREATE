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


process seqtk_subset
{
    publishDir "${params.outdir}/filtered", mode: 'copy', overwrite: params.force, enable: params.keep
    container "alexdhill/create:seqtk"
    if (params.manage_resources)
    {
        cpus 8
        memory '16.GB'
    }
    input:
        tuple(
            val(sample),
            val(nreads),
            path(controls),
            path(read)
        )
    output:
        tuple(
            val("${sample}"),
            env(NREADS),
            path("${sample}_filtered.fq.gz")
        )
    shell:
        '''
            mkfifo reads read_names control_names
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Removing DCS Sequence(s)"
                echo "Read: !{read}"
                echo "Controls: !{controls}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            cat !{controls} \
            | sort \
            > control_names \
            & zcat !{read} \
            | sed -n '1~4p' \
            | cut -d' ' -f1 \
            | sed 's/@//' \
            | sort \
            > read_names \
            & comm -23 read_names control_names \
            > non_control_reads.txt \

            gzip -cd !{read} > reads &
            seqtk subseq reads non_control_reads.txt \
            | gzip \
            > !{sample}_filtered.fq.gz

            NREADS="$(cat non_control_reads.txt | wc -l)"
        '''
}