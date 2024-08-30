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


include { alevin_align } from "../../../modules/salmon/alevin_align/alevin_align.nf"
include { alevin_collate } from "../../../modules/salmon/alevin_collate/alevin_collate.nf"
include { alevin_quantify } from "../../../modules/salmon/alevin_quantify/alevin_quantify.nf"
include { compile_quantifications } from "../../../modules/R/compile_quantifications/compile_quantifications.nf"
include { run_analysis } from "../../../modules/R/run_analysis/run_analysis.nf"

workflow SINGLE_CELL
{
    take:
        reads
    main:
        if (params.library=="single_cell")
        {
            reference = Channel.fromPath(params.ref)
            barcodes = Channel.fromPath(params.barcodes)

            reads
            | combine(reference)
            | alevin_align
            | combine(barcodes)
            | alevin_collate
            | combine(reference)
            | alevin_quantify
            | collect
            | map{quants -> [quants]}
            | combine(reference)
            | compile_quantifications
        }
}