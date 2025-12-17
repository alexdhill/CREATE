
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


process minimap2_align_genome
{
    publishDir "${params.outdir}/align", mode: 'copy', overwrite: params.force, enabled: params.keep
    container 'alexdhill/create:minimap2-2.26'
    conda projectDir+'/bin/conda/modules/minimap2.yaml'
    if (params.manage_resources)
    {
        cpus 8
        memory '32.GB'
    }
    input:
        tuple(
            val(sample),
            val(nreads),
            path(read),
            path(reference),
            path(parameters)
        )
    output:
        tuple(
            val("${sample}"),
            val("${nreads}"),
            path("${sample}.coord_sorted.bam")
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Aligning to Minimap2 Index"
                echo "Sample: !{sample}"
                echo "Read: !{read} (!{nreads} reads)"
                echo "Reference: !{reference}"
                echo "User parameters: !{parameters}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            params="--eqx -N 10000"
            if [[ "!{parameters}" != "NULL" ]]; then
                params="$(jq '.minimap2 | to_entries | .[] | "\\(.key)=\\(.value)"' !{parameters} | xargs | sed 's/=true//g')"
            fi

            minimap2 -ax splice:sr \
                ${params} \
                -t !{task.cpus} \
                !{reference}/*genome.fa.gz \
                !{read} \
            > !{sample}_raw.bam
            samtools sort !{sample}_raw.bam \
            | samtools view -b - \
            > !{sample}.coord_sorted.bam
            '''
}