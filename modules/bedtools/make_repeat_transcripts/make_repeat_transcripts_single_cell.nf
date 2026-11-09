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
 * Repeat transcripts for the single-cell (splintr) index.
 *
 * The splintr index does not use the bulk repeat FASTA: make_splintr_transcripts
 * extracts every target from the genome by its GTF coordinates
 * (eisaR + extractTranscriptSeqs). The flank therefore has to be applied to the
 * repeat ANNOTATION that feeds it. This process writes a copy of the repeat GTF
 * whose gene/transcript/exon records are extended by --flank_single_cell nt
 * (clipped at contig ends). transcript_id/gene_id are unchanged, so the
 * _range=..._strand= locus names, tx2g/tx3g maps and host tags keep pointing at
 * the element itself.
 *
 * --flank_scope_single_cell:
 *   nonexonic (default)  flank only loci that do not overlap a GENCODE exon
 *                        (either strand). In 3' single-cell data an exonic
 *                        (e.g. 3' UTR Alu) flank takes host-gene UMIs.
 *   all                  flank every locus
 *
 * Downstream, the intron mask is computed from these (flanked) repeat ranges,
 * so flank sequence is removed from intron targets as well and is not
 * offered twice.
 *
 * The flanked GTF is used ONLY to build the splintr transcripts; the published
 * complete annotation, linked txome and transcript map keep element
 * coordinates.
 */
process make_repeat_transcripts_single_cell
{
    publishDir "${params.outdir}/", mode: 'copy', enabled: params.keep, overwrite: params.force
    container 'alexdhill/create:bedtools-2.31.0'
    conda projectDir+'/bin/conda/modules/bedtools.yaml'
    if (params.manage_resources)
    {
        cpus 4
        memory '16.GB'
    }
    input:
        tuple(
            path(repeat_gtf),
            path(annotation), // GENCODE GTF, or assets/NULL for a repeat-only build
            path(reference)
        )
    output:
        path("*_single_cell_repeat_annotation.gtf.gz")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Creating single-cell repeat annotation..."
                echo "Repeats: !{repeat_gtf}"
                echo "Flank: !{params.flank_single_cell} nt (!{params.flank_scope_single_cell})"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi
            set -o pipefail

            version='!{params.version}'
            if [[ "!{params.genome}" == "T2T" ]]; then
                version='2'
            fi

            # Repeat GTFs are 1-based (start = BED start + 1). Guard only: bedtools
            # rejects GTF start < 1, which a 0-based GTF would produce at T2T's
            # chromosome-initial telomeric (CTAACC)n loci. No-op on correct input.
            zcat -f !{repeat_gtf} \
            | awk -F'\\t' -v OFS='\\t' 'NF>=9 { if ($4 < 1) $4 = 1; print }' \
            > repeats.gtf
            flank=!{params.flank_single_cell}

            if [[ "${flank}" -gt 0 ]]; then
                pigz -cdp !{task.cpus} !{reference} \
                | awk '/^>/ {if (seqlen) print chr "\\t" seqlen; chr=substr($1,2); seqlen=0; next} {seqlen += length($0)} END {if (seqlen) print chr "\\t" seqlen}' \
                > chrom_sizes.tsv

                if [[ "!{params.flank_scope_single_cell}" == "nonexonic" && "$(basename !{annotation})" != "NULL" ]]; then
                    zcat -f !{annotation} \
                    | awk -F'\\t' -v OFS='\\t' '$3=="exon" {print $1, $4-1, $5}' \
                    > exons.bed
                    # gene/transcript/exon records of a locus share coordinates,
                    # so all three land on the same side of this split
                    bedtools intersect -u -a repeats.gtf -b exons.bed > exonic.gtf
                    bedtools intersect -v -a repeats.gtf -b exons.bed \
                    | bedtools slop -i stdin -g chrom_sizes.tsv -b ${flank} \
                    | cat - exonic.gtf \
                    > flanked.gtf
                    echo "Flanked $(( ($(wc -l < repeats.gtf) - $(wc -l < exonic.gtf)) / 3 )) non-exonic loci; $(( $(wc -l < exonic.gtf) / 3 )) exonic loci unflanked"
                else
                    bedtools slop -i repeats.gtf -g chrom_sizes.tsv -b ${flank} > flanked.gtf
                fi
            else
                cp repeats.gtf flanked.gtf
            fi

            [[ "$(wc -l < flanked.gtf)" -eq "$(wc -l < repeats.gtf)" ]] \
                || { echo "ERR: flanking changed the number of repeat records" 1>&2; exit 1; }

            gzip -c flanked.gtf > !{params.genome}v${version}_single_cell_repeat_annotation.gtf.gz
        '''
}
