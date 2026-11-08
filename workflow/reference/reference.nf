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
import groovy.util.logging.Slf4j

include { download_reference } from '../../modules/bash/download_reference/download_reference.nf'
include { download_gencode_transcripts } from '../../modules/bash/download_gencode_transcripts/download_gencode_transcripts.nf'
include { download_gencode_annotation } from '../../modules/bash/download_gencode_annotation/download_gencode_annotation.nf'
include { correct_isoquant } from '../../modules/bash/correct_isoquant/correct_isoquant.nf'
include { generate_transcriptome } from '../../modules/bash/generate_transcriptome/generate_transcriptome.nf'
include { download_repeat_regions } from '../../modules/bash/download_repeat_regions/download_repeat_regions.nf'
include { make_repeat_transcripts_short } from '../../modules/bedtools/make_repeat_transcripts/make_repeat_transcripts_short.nf'
include { make_repeat_transcripts_long } from '../../modules/bedtools/make_repeat_transcripts/make_repeat_transcripts_long.nf'
include { make_repeat_transcripts_single_cell } from '../../modules/bedtools/make_repeat_transcripts/make_repeat_transcripts_single_cell.nf'
include { download_repeat_annotation } from '../../modules/R/download_repeat_annotation/download_repeat_annotation.nf'
include { make_complete_transcripts } from '../../modules/bash/make_complete_transcripts/make_complete_transcripts.nf'
include { make_complete_annotation } from '../../modules/bash/make_complete_annotation/make_complete_annotation.nf'
include { make_complete_transcripts as make_complete_transcripts_long } from '../../modules/bash/make_complete_transcripts/make_complete_transcripts.nf'
include { make_complete_annotation as make_complete_annotation_single_cell } from '../../modules/bash/make_complete_annotation/make_complete_annotation.nf'
include { make_transcript_map } from '../../modules/R/make_transcript_map/make_transcript_map.nf'
include { link_transcriptome } from '../../modules/R/link_transcriptome/link_transcriptome.nf'
include { make_database } from '../../modules/R/make_database/make_database.nf'

include { SHORT } from '../../subworkflow/reference/short/short.nf'
include { LONG } from '../../subworkflow/reference/long/long.nf'
include { SINGLE_CELL } from '../../subworkflow/reference/single_cell/single_cell.nf'
include { DISCOVER } from '../../subworkflow/reference/discover/discover.nf'

workflow REFERENCE
{
    noent = Channel.fromPath(projectDir+'/assets/NULL')
    
    log.info("Starting main logic...")
    // Get reference genome
    download_reference()
    | set{reference}

    // Get gene annotation and transcriptome
    if (!params.reonly) {
        if (params.isoquant) {
            log.info("Using provided isoquant annotation")
            Channel.fromPath(params.isoquant)
            | correct_isoquant
            | set{annotation}

            annotation
            | combine(reference)
            | generate_transcriptome
            | set{transcriptome}
        } else {
            log.info("Getting GENCODE annotation and transcripts")
            download_gencode_annotation()
            | combine(reference)
            | download_gencode_transcripts
            | set{transcriptome}

            download_gencode_annotation.out
            | set{annotation}
        }
    }

    reference
    | DISCOVER

    if (params.txonly) {
        log.info("Making repeat-naive reference directory...")
        annotation
        | combine(transcriptome)
        | combine(reference)
        | link_transcriptome & SHORT & LONG & SINGLE_CELL

        annotation
        | combine(noent)
        | make_transcript_map
    } else if (params.reonly) {
        log.info("Making repeat-only reference directory...")
        // Get RepeatMasker track
        download_repeat_regions()
        | set{repeatmasker}

        reference
        | combine(repeatmasker) // [genome.fa, repeats.bed]
        | download_repeat_annotation
        | set{repeat_gtf}

        // Repeat sequences are built per platform (see make_repeat_transcripts_*):
        // the flank each index wants differs. There is no GENCODE annotation here,
        // so a 'nonexonic' scope flanks every locus.
        reference
        | combine(repeatmasker)
        | combine(noent) // [genome.fa, repeats.bed, NULL]
        | make_repeat_transcripts_short
        | set{repeat_fa_short}

        repeat_gtf
        | combine(repeat_fa_short)
        | combine(reference) // [repeats.gtf, repeats_short.fa, genome.fa]
        | set{complete_short}

        complete_short
        | link_transcriptome & SHORT

        if (index_requested('long')) {
            reference
            | combine(repeatmasker)
            | combine(noent)
            | make_repeat_transcripts_long
            | set{repeat_fa_long}

            repeat_gtf
            | combine(repeat_fa_long)
            | combine(reference) // [repeats.gtf, repeats_long.fa, genome.fa]
            | LONG
        }

        if (index_requested('single_cell')) {
            repeat_gtf
            | combine(noent)
            | combine(reference) // [repeats.gtf, NULL, genome.fa]
            | make_repeat_transcripts_single_cell // flanked repeats.gtf
            | combine(repeat_fa_short)
            | combine(reference) // [repeats_sc.gtf, (unused) repeats.fa, genome.fa]
            | SINGLE_CELL
        }

        noent
        | combine(repeat_gtf)
        | make_transcript_map
    } else {
        log.info("Making repeat-aware reference directories...")
        // Get RepeatMasker track
        download_repeat_regions()
        | set{repeatmasker}

        reference
        | combine(repeatmasker) // [genome.fa, repeats.bed]
        | download_repeat_annotation
        | set{repeat_gtf}

        // Short read: flanked repeat FASTA (--flank_short / --flank_scope_short).
        // This complete FASTA is also the one the linked txome is digested from.
        reference
        | combine(repeatmasker)
        | combine(annotation) // [genome.fa, repeats.bed, gencode.gtf]
        | make_repeat_transcripts_short
        | set{repeat_fa_short}

        repeat_gtf
        | combine(repeat_fa_short)
        | combine(annotation)
        | combine(transcriptome) // [repeats.gtf, repeats_short.fa, gencode.gtf, gencode.fa]
        | (make_complete_annotation & make_complete_transcripts)

        make_complete_annotation.out
        | combine(make_complete_transcripts.out)
        | combine(reference) // [complete.gtf, complete_short.fa, genome.fa]
        | set{complete_short}

        complete_short
        | link_transcriptome & SHORT

        // Long read: its own repeat FASTA (--flank_long, default 0).
        if (index_requested('long')) {
            reference
            | combine(repeatmasker)
            | combine(annotation)
            | make_repeat_transcripts_long
            | set{repeat_fa_long}

            repeat_gtf
            | combine(repeat_fa_long)
            | combine(annotation)
            | combine(transcriptome)
            | make_complete_transcripts_long

            make_complete_annotation.out
            | combine(make_complete_transcripts_long.out)
            | combine(reference) // [complete.gtf, complete_long.fa, genome.fa]
            | LONG
        }

        // Single cell: splintr targets are cut from the genome by GTF coordinates,
        // so the flank (--flank_single_cell / --flank_scope_single_cell) is applied
        // to a private copy of the repeat GTF used only for this index.
        if (index_requested('single_cell')) {
            repeat_gtf
            | combine(annotation)
            | combine(reference) // [repeats.gtf, gencode.gtf, genome.fa]
            | make_repeat_transcripts_single_cell
            | combine(repeat_fa_short)
            | combine(annotation)
            | combine(transcriptome) // [repeats_sc.gtf, (unused) fa, gencode.gtf, gencode.fa]
            | make_complete_annotation_single_cell
            | combine(make_complete_transcripts.out)
            | combine(reference) // [complete_sc.gtf, (unused) fa, genome.fa]
            | SINGLE_CELL
        }

        annotation
        | combine(repeat_gtf)
        | make_transcript_map
    }
}

// The long and single-cell repeat builds are skipped unless their index is asked
// for (DISCOVER needs the long index).
def index_requested(name) {
    def wanted = params.index.toString().split(',')
    return wanted.contains(name) || (name == 'long' && wanted.contains('discover'))
}
