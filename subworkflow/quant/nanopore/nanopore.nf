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


include { count_reads_np } from "../../../modules/bash/count_reads/count_reads_np.nf"
include { minimap2_align_dcs } from "../../../modules/minimap2/minimap2_align/minimap2_align_dcs.nf"
include { trim_reads_nanopore } from "../../../modules/nanoplot/trim_reads/trim_reads_nanopore.nf"
include { minimap2_align } from "../../../modules/minimap2/minimap2_align/minimap2_align.nf"
include { salmon_quant_nanopore } from "../../../modules/salmon/salmon_quant/salmon_quant_nanopore.nf"
include { compile_quantifications } from "../../../modules/R/compile_quantifications/compile_quantifications.nf"

workflow NANOPORE
{
    take:
        is_acc
    main:
        if (params.library=="nanopore")
        {
            reference = Channel.fromPath(params.ref)
            dcs = Channel.fromPath(params.dcs)
            parameters = Channel.fromPath(params.parameters)
            metadata = Channel.fromPath(params.metadata)
            if (is_acc)
            {
                log.error("Nanopore data cannot be hosted from SRA. Exiting")
                error(2)
            } else
            {
                reads = Channel.fromPath(params.samples+"/"+params.pattern)
                .map{sample -> [sample.name.split(/\.f(ast)?q(\.gz)?/)[0], sample]}
            }

            reads
            | count_reads_np
            | concat(reads)
            | groupTuple()
            | map{ res -> [res[0], res[1][0], res[1][1]] }
            | combine(dcs)
            | minimap2_align_dcs
            | trim_reads_nanopore
            | combine(reference)
            | combine(parameters)
            | minimap2_align
            | combine(reference)
            | combine(parameters)
            | salmon_quant_nanopore
            | collect
            | map{quants -> [quants]}
            | combine(reference)
            | combine(metadata)
            | compile_quantifications
        }
}