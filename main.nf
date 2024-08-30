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
 * Imports
 */
import groovy.util.logging.Slf4j
// ERROR -> WARN -> DEBUG -> INFO


/*
 * Print final parameters
 */
def print_val(k,v)
{
    type=v.getClass().getSimpleName()
    key = k
    while (key.length()<12) key += " "
    switch(type)
    {
        case "Boolean":
            log.info("    $key => \u001B[32m$v\u001B[0m")
            return
        case "String":
            value = v
            value = value.replace(',', '\u001B[0m,\u001B[33m')
            value = value.replace('*', '\u001B[35m*\u001B[33m')
            log.info("    $key => \u001B[33m$value\u001B[0m")
            return
        case "Integer":
            log.info("    $key => \u001B[34m$v\u001B[0m")
            return
        case "Float":
            log.info("    $key => \u001B[34m$v\u001B[0m")
            return
        default:
            log.info("    $key => $v")
            return
    }
}
log.info("Run parameters:")
params.each{k, v -> ['',false,0].contains(v)?:print_val(k,v)}
log.info(" ")
params.manage_resources = (params.limits || params.exec!='local')

/*
 * Handle post-print params
 */
if (params.get('reference'))
{
    if (![null,false].contains(params.get('isoquant'))) params.version = '00'
    else if (params.get('genome')=='T2T') params.version = '2'
}


/*
 * Make CREATE workflow selector
 */
include { REFERENCE } from './workflow/reference/reference.nf'
include { QUANT } from './workflow/quant/quant.nf'

workflow CREATE
{
    if (params.get('quant')!=null && params.get('quant')==true) QUANT()
    else REFERENCE()
}


/*
 * Run CREATE
 */
workflow { CREATE() }