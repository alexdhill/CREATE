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
 

process download_repeat_annotation
{
    publishDir "${params.outdir}/", mode: 'copy', enabled: params.keep, overwrite: params.force
    if (params.manage_resources)
    {
        cpus 1
        memory '32.GB'
    }
    input:
        tuple(
            path(genome),
            path(regions)
        )
    output:
        path("*.gtf.gz")
    shell:
        if (params.genome=="T2T")
        {
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Generating T2T repeat annotation..."
                echo "Repeats: !{regions}"
                echo "Reference: !{genome}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set +ex
            fi

            cat !{regions} \
            | awk -F"\\t" '{print \
$1"\\tT2T_rmsk\\tgene\\t"$2"\\t"$3"\\t"$5"\\t"$6"\\t.\\tgene_id \\""$4"\\"; gene_name \\""$4"\\"; transcript_id \\""$4"_range="$1":"$2"-"$3"_strand="$6"\\"; gene_biotype \\""$8","$7"\\"; \\n"\
$1"\\tT2T_rmsk\\ttranscript\\t"$2"\\t"$3"\\t"$5"\\t"$6"\\t.\\tgene_id \\""$4"\\"; gene_name \\""$4"\\"; transcript_id \\""$4"_range="$1":"$2"-"$3"_strand="$6"\\"; gene_biotype \\""$8","$7"\\"; \\n"\
$1"\\tT2T_rmsk\\texon\\t"$2"\\t"$3"\\t"$5"\\t"$6"\\t0\\tgene_id \\""$4"\\"; gene_name \\""$4"\\"; transcript_id \\""$4"_range="$1":"$2"-"$3"_strand="$6"\\";  gene_biotype \\""$8","$7"\\"; exon_numer \\"1\\"; "}' \
            | gzip --best \
            > T2Tv2_repeat_annotation.gtf.gz
            ## GTF: seq source feature start stop score strand frame group
            ## BED: chr start stop locus score strand super family ?? ?? 
        '''
        }
        else
        {
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Downloading HG38 repeat annotation..."
            fi
            verbose=""
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
                verbose="-v"
            fi

            bash !{projectDir}/bin/sql/repeatmasker_info.mysql > repeat_info.tsv
            bash !{projectDir}/bin/sql/repeatmasker_annotation.mysql > raw_rmsk.gtf
            for chr in `seq 1 22` X Y; do
                cat raw_rmsk.gtf | grep  "^chr${chr}\t" | sed -re"s/\t$//"
            done \
            > filtered.gtf

            sed -re's/^.+gene_id \"//' -e's/\".+$//' filtered.gtf \
            > rep_list.txt

            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Reformatting annotation..."
            fi
            Rscript !{projectDir}/bin/R/make_repeat_annotation.R \
                --gtf filtered.gtf \
                --info repeat_info.tsv \
                --list rep_list.txt \
                --out filtered_biotyped.gtf

            sed 's/""/"/g' filtered_biotyped.gtf \
            | awk 'match($0, /gene_id "([^"]+)"/, a) { print $0 " gene_name \"" a[1] "\"; " }' \
            > HG38v!{params.version}_repeat_annotation.gtf

            if [[ `wc -l HG38v!{params.version}_repeat_annotation.gtf` == '0' ]]; then
                echo "Error formatting GTF"
                exit 1
            fi
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Compressing annotation..."
            fi
            gzip HG38v!{params.version}_repeat_annotation.gtf
        '''
        }
}
