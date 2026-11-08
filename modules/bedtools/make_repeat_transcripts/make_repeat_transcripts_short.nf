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

/*
 * Repeat transcripts for the short index.
 *
 * Each repeat locus is extended by --flank_short nt on both sides
 * (clipped at contig ends) before its sequence is extracted, so reads whose
 * fragment runs past the element boundary still have a full-length target.
 * Locus names are kept unchanged, so the transcript map is shared by all
 * indexes.
 *
 * --flank_scope_short chooses which loci are flanked:
 *   all        every locus
 *   nonexonic  only loci that do not overlap a GENCODE exon (either strand);
 *              exonic loci keep their element coordinates so host mRNA
 *              fragments are not offered a flanked repeat target.
 *
 * The flank is a sequence-only change: annotation coordinates, locus names and
 * the transcript map are untouched.
 */
process make_repeat_transcripts_short
{
    publishDir "${params.outdir}/", mode: 'copy', enabled: params.keep, overwrite: params.force
    container 'alexdhill/create:bedtools-2.31.0'
    conda projectDir+'/bin/conda/modules/bedtools.yaml'
    if (params.manage_resources)
    {
        cpus 8
        memory '16.GB'
    }
    input:
        tuple(
            path(reference),
            path(regions),
            path(annotation) // GENCODE GTF, or assets/NULL for a repeat-only build
        )
    output:
        path("*_short_repeat_transcripts.fa.gz")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Creating short repeat transcripts..."
                echo "Repeats: !{regions}"
                echo "Reference: !{reference}"
                echo "Flank: !{params.flank_short} nt (!{params.flank_scope_short})"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi
            set -o pipefail

            if [[ "!{params.genome}" == "T2T" ]]; then
                version='2'
                cat !{regions} \
                | awk '{print $1"\\t"$2"\\t"$3"\\t"$4"_range="$1":"$2"-"$3"_strand="$6"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10}' \
                > regions.bed
            else
                version='!{params.version}'
                cp !{regions} regions.bed
            fi

            pigz -cdp !{task.cpus} !{reference} > genome
            awk '/^>/ {if (seqlen) print chr "\\t" seqlen; chr=substr($1,2); seqlen=0; next} {seqlen += length($0)} END {if (seqlen) print chr "\\t" seqlen}' genome \
            > chrom_sizes.tsv

            flank=!{params.flank_short}
            if [[ "${flank}" -gt 0 && "!{params.flank_scope_short}" == "nonexonic" && "$(basename !{annotation})" != "NULL" ]]; then
                zcat -f !{annotation} \
                | awk -F'\\t' -v OFS='\\t' '$3=="exon" {print $1, $4-1, $5}' \
                > exons.bed
                bedtools intersect -u -a regions.bed -b exons.bed > exonic.bed
                bedtools intersect -v -a regions.bed -b exons.bed \
                | bedtools slop -i stdin -g chrom_sizes.tsv -b ${flank} \
                | cat - exonic.bed \
                > flanked.bed
                echo "Flanked $(( $(wc -l < regions.bed) - $(wc -l < exonic.bed) )) non-exonic loci; $(wc -l < exonic.bed) exonic loci unflanked"
            elif [[ "${flank}" -gt 0 ]]; then
                bedtools slop -i regions.bed -g chrom_sizes.tsv -b ${flank} > flanked.bed
            else
                cp regions.bed flanked.bed
            fi

            # every locus must survive, exactly once
            [[ "$(wc -l < flanked.bed)" -eq "$(wc -l < regions.bed)" ]] \
                || { echo "ERR: flanking changed the number of repeat loci" 1>&2; exit 1; }

            bedtools getfasta -nameOnly -fi genome -bed flanked.bed \
            | pigz -cp !{task.cpus} \
            > !{params.genome}v${version}_short_repeat_transcripts.fa.gz
        '''
}
