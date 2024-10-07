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
 

include { salmon_index_novel } from "../../../modules/salmon/salmon_index/salmon_index_novel.nf"
include { minimap2_index_novel } from "../../../modules/minimap2/minimap2_index/minimap2_index_novel.nf"

workflow NOVEL_INDEX
{
    take:
        isoforms

    main:
        isoforms
        | map{isoforms -> [isoforms[0], isoforms[4]]}
        | salmon_index_novel & minimap2_index_novel

    emit:
        short_index=salmon_index_novel.out
        long_index=minimap2_index_novel.out

}