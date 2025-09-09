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


process oarfish_quant
{
    publishDir "${params.outdir}/quant/", mode: 'copy', enabled: params.keep, overwrite: params.force
    container 'alexdhill/create:oarfish-0.9.0'
    conda projectDir+'/bin/conda/modules/oarfish.yaml'
    if (params.manage_resources)
    {
        cpus 8
        memory '32.GB'
        time '6h'
    }
    input:
        tuple(
            val(sample),
            val(nreads),
            path(bam),
            path(reference),
            path(parameters)
        )
    output:
        path("${sample}/")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Quantifying nanopore reads"
                echo "Sample: !{sample} (!{nreads} reads)"
                echo "BAM: !{bam}"
                echo "Reference: !{reference}"
                echo "User parameters: $(jq '.salmon' !{parameters})"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            params="--filter-group no-filters --model-coverage"
            if [[ "${parameters}" != "NULL" ]]; then
                params="$(jq '.oarfish | to_entries | .[] | "--\\(.key)=\\(.value)"' flags.json | xargs | sed 's/=true//g')"
            fi

            oarfish quant \
                ${params} \
                -a !{bam} \
                -j !{task.cpus} \
                -o !{sample}/!{sample}

            if [[ ! -e !{sample}/!{sample}.ambig_info.tsv ]]; then
                echo "\033[1;31mERR: Oarfish quantification failed\033[0m" 1>&2
                exit 1
            fi

            grep 'SeqHash' !{reference}/*complete_digest/info.json \
            | awk '{print $2}' \
            | sed -r 's/[\",]//g' \
            | xargs -I'%' \
                sed -ri \
                    's/\"sha256_seqs\": \"[a-zA-Z0-9]+\"/\"sha256_seqs\": \"%\"/' \
                    !{sample}/!{sample}.meta_info.json

            grep 'NameHash' !{reference}/*complete_digest/info.json \
            | awk '{print $2}' \
            | sed -r 's/[\",]//g' \
            | xargs -I'%' \
                sed -ri \
                    's/\"sha256_names\": null/\"sha256_names\": \"%\"/' \
                    !{sample}/!{sample}.meta_info.json
        '''
}