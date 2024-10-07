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

import java.nio.file.FileSystems
import java.io.File
import java.util.Arrays

include { count_reads_pe } from "../../modules/bash/count_reads/count_reads_pe.nf"
include { count_reads_np } from "../../modules/bash/count_reads/count_reads_np.nf"
include { minimap2_align_dcs } from "../../modules/minimap2/minimap2_align/minimap2_align_dcs.nf"
include { seqtk_subset } from "../../modules/seqtk/seqtk_subset/seqtk_subset.nf"
include { trim_reads_nanopore } from "../../modules/nanoplot/trim_reads/trim_reads_nanopore.nf"
include { trim_reads_paired } from "../../modules/trim-galore/trim_reads/trim_reads_paired.nf"
include { star_align_genome } from "../../modules/star/star_align/star_align_genome.nf"
include { flair_junctions } from "../../modules/flair/flair_junctions/flair_junctions.nf"
include { flair_align } from "../../modules/flair/flair_align/flair_align.nf"
include { flair_correct } from "../../modules/flair/flair_correct/flair_correct.nf"
include { flair_collapse } from "../../modules/flair/flair_collapse/flair_collapse.nf"
include { make_novel_reference } from "../../modules/bash/make_novel_reference/make_novel_reference.nf"
include { NOVEL_INDEX } from "../../subworkflow/discover/novel_index/novel_index.nf"
include { minimap2_align } from "../../modules/minimap2/minimap2_align/minimap2_align.nf"
include { salmon_quant_novel } from "../../modules/salmon/salmon_quant/salmon_quant_novel.nf"
include { compile_quantifications } from "../../modules/R/compile_quantifications/compile_quantifications.nf"

workflow DISCOVER
{
    long_dir = (params.long_reads.lastIndexOf("/")+1==params.long_reads.length())?params.long_reads:params.long_reads+"/"
    long_glob = FileSystems.getDefault().getPathMatcher("glob:${long_dir}${params.lr_pattern}")
    long_files = Arrays.stream(new File(long_dir).listFiles())
        .filter(f -> long_glob.matches(f.toPath()))
        .toArray();
    Arrays.sort(long_files);
    log.info("Listing long reads: ${long_dir}${params.lr_pattern}")
    Arrays.copyOfRange(long_files, 0, Math.min(long_files.length, 10)).eachWithIndex{f,i -> {
        if (i%2<1) log.info(">   \u001B[36m${f.toString()}\u001B[0m")
        else log.info(">   ${f.toString()}")
    }}
    if (long_files.size() > 10) log.info(">   . . . and \u001B[32m${long_files.size()-10}\u001B[0m more")
    log.info(" ")
    
    paired_dir = (params.paired_reads.lastIndexOf("/")+1==params.paired_reads.length())?params.paired_reads:params.paired_reads+"/"
    paired_glob = FileSystems.getDefault().getPathMatcher("glob:${paired_dir}${params.pr_pattern}")
    paired_files = Arrays.stream(new File(paired_dir).listFiles())
        .filter(f -> paired_glob.matches(f.toPath()))
        .toArray()
    Arrays.sort(paired_files)
    log.info("Listing paired reads: ${paired_dir}${params.pr_pattern}")
    Arrays.copyOfRange(paired_files, 0, Math.min(paired_files.length, 10)).eachWithIndex{f,i -> {
        if (i%4<2) log.info(">   \u001B[36m${f.toString()}\u001B[0m")
        else log.info(">   ${f.toString()}")
    }}
    if (paired_files.size() > 10) log.info(">   . . . and \u001B[32m${paired_files.size()-10}\u001B[0m more")
    log.info(" ")

    if (paired_files.length/2 != long_files.length)
    {
        log.warn("Number of long reads and paired read samples do not match")
    }

    prefix_list = (new File(params.prefix) as String[]).sort(false)
    log.info("Found ${prefix_list.length} prefixes")
    for (prefix in prefix_list)
    {
        if (long_files.find{f -> f.name.startsWith(prefix)} == null)
        {
            log.error("No long read sample found for prefix: ${prefix}")
            throw new Exception("Prefix ${prefix} missing samples")
        }
        if (paired_files.find{f -> f.name.startsWith(prefix)} == null)
        {
            log.error("No paired read sample found for prefix: ${prefix}")
            throw new Exception("Prefix ${prefix} missing samples")
        }
    }

    reference = Channel.fromPath(params.ref)
    dcs = Channel.fromPath(params.dcs)
    prefixes = Channel.from(prefix_list)

    paired_channel = Channel.fromFilePairs(paired_dir+params.pr_pattern)
        .map{sample -> [sample[0], sample[1][0], sample[1][1]]}
        .toSortedList{a,b -> a[0] <=> b[0]}
        .flatMap{it}
        .merge(prefixes)
        .map{sample -> [sample[3], sample[1], sample[2]]}

    count_reads_pe(paired_channel)
    | trim_reads_paired
    | combine(reference)
    | star_align_genome
    | map{sample -> [sample[0], sample[1], sample[2]]}
    | combine(reference)
    | flair_junctions

    long_channel = Channel.fromPath(long_dir+params.lr_pattern)
        .map{file -> [file.name, file]}
        .toSortedList{a,b -> a[0] <=> b[0]}
        .flatMap{it}
        .merge(prefixes)
        .map{sample -> [sample[2], sample[1]]}

    count_reads_np(long_channel)
    | combine(dcs)
    | minimap2_align_dcs
    | join(long_channel)
    | seqtk_subset
    | trim_reads_nanopore
    | combine(reference)
    | flair_align
    | join(flair_junctions.out)
    | combine(reference)
    | flair_correct
    | map{res -> res[1]}
    | collect
    | map{beds -> [beds]}
    | combine(
        trim_reads_nanopore.out
        | map{res -> res[2]}
        | collect
        | map{reads -> [reads]}
    )
    | combine(reference)
    | flair_collapse
    | combine(reference)
    | make_novel_reference & NOVEL_INDEX

    trim_reads_paired.out
    | combine(NOVEL_INDEX.out.short_index)
    | salmon_quant_novel
    | collect
    | map{quants -> [quants]}
    | combine(reference)
    | compile_quantifications
}
