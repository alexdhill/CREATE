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
 * Handle parameters
 */

params.help = false
if (params.help)
{
    println '''Usage: nextflow run alexdhill/CREATE [--quant/--reference/--discover] [options]
    reference         Build a reference for CREATE
    quant             Quantify RNA-seq reads with CREATE
    discover          Discover novel isoforms [BETA]

options:
    Globals:
        --help        Print this help message
        --container   Containerization method [docker(default),conda] TODO!
        --outdir      Output directory [REQUIRED]
        --log         Log level [OFF(default),INFO,DEBUG]
        --exec        Executor mode [local(default),slurm]
        --force       Overwrite existing files [default=false]
        --keep        Keep intermediate files [default=false]
        --limits      Add sample-based resources limits when running locally [default=false]
        --clean       Clean up working directory (breaks `-resume`!) [default=false]

    --reference
        --genome      Genome assembly [HG38,T2T] [REQUIRED]
        --version     Gencode version [HG38 only] [default=39]
        --index       Desired index [short,long,single_cell] [default=short]
        --isoquant    Novel isoform annotation [BETA]

    --quant
        --samples     Input samples directory [REQUIRED]
        --ref         CREATE reference directory [REQUIRED]
        --library     Library type [paired_end,nanopore,single_cell] [REQUIRED]
                      NOTE: the reference provided must contain a compatible index
        --pattern     File pattern for input samples
        --metadata    An unnamed, two-column CSV file of samples,condition

    --discover
        --ref          CREATE reference directory [REQUIRED]
        --dcs          The DNA/RNA control sequence fasta file [REQUIRED]
        --prefixes     A list of file prefixes for the paired long-short reads [REQUIRED]
        --long_reads   The directory of the long-read sequences [REQUIRED]
        --paired_reads The directory of the paired-end sequences [REQUIRED]
        --lr_pattern   The pattern for the long-read files [default='*.fastq.gz']
        --pr_pattern   The pattern for the paired-end files [default='*_R{1,2}_*.fastq.gz']
        --dump         Location to save the 'novel' CREATE reference [default=false]

    Subworkflow-specific options:
        [--exec local]
            --threads     Maximum number of threads to allocate TODO!
            --memory      Maximum memory to allocate TODO!
        [--exec slurm]
            --njobs       Maximum number of jobs to run concurrently [default=15]
            --scratch     Scratch directory for temporary files [default=false]
            --account     The user account to submit jobs
        [--quant --library single_cell]
            --barcodes    A list (txt) of all barcodes used in the experiment [REQUIRED]
            --chemistry   The single-cell chemistry used [dropseq,chromium,chromiumV3(default)]
        [--quant --library nanopore]
            --dcs         File with the DCS/RCS sequence used (for removal) [REQUIRED]
'''
    System.exit 0
}

params.outdir = ''
params.log = 'OFF'
params.keep = false
params.limits = false
params.threads = 0
params.memory = 0
params.force = false
params.quant = false
params.reference = false
if (params.reference)
{
    params.genome = ''
    params.isoquant = false
    params.index = ''
    if (params.genome=='HG38' || params.genome=="MM39")
    {
        params.version = 39
        if (params.version < 1)
        {
            println("--version: ${params.version}")
            throw new IllegalArgumentException('Gencode version must be at least 1')
            if ((params.genome[0]=="M") && !(params.version[0]=="M"))
            {
                println("--genome: ${params.genome}")
                println("--version: ${params.version}")
                throw new IllegalArgumentException('Gencode version MUST start with `M` when building mouse reference')
            }
        }
    } else if (params.genome=="T2T")
    {
        params.version = -1
        if (params.version!=-1)
        {
            println("--version: ${params.version}");
            throw new IllegalArgumentException('Gencode version cannot be specified for T2T genome')
        }
    }
    else if (params.genome!="T2T")
    {
        println("--genome: ${params.genome}");
        throw new IllegalArgumentException('Genome must be either "HG38" or "T2T"')
    }
    index = new ArrayList(Arrays.asList(params.index.split(',')))
    index.retainAll(['short','long','single_cell','discover'])
    if (index.size()!=params.index.split(',').size())
    {
        println("--index: ${params.index}");
        throw new IllegalArgumentException('index must be one or more of: short, long, single_cell')
    }
}
else if (params.quant)
{
    params.samples = ''
    params.ref = ''
    params.library = 'paired_end'
    params.metadata = ''
    if (['paired_end','single_cell'].contains(params.library))
    {
        params.pattern = '*_R{1,2}_*.fastq.gz'
        if (params.library=="single_cell")
        {
            params.barcodes = ''
            params.chemistry = "chromiumV3"
            if (params.barcodes=='')
            {
                println("--barcodes: ${params.barcodes}");
                throw new IllegalArgumentException('Barcode list must be specified')
            }
            if (!["chromiumV3","chromium","dropseq"].contains(params.chemistry))
            {
                println("--chemistry: ${params.chemistry}");
                throw new IllegalArgumentException('Chemistry must be one of: chromiumV3, chromium, dropseq')
            }
        }
    }
    else if (params.library=='nanopore')
    {
        params.pattern = '*.fastq.gz'
        params.dcs = ""
        params.discovery = false
        if (params.dcs=='')
        {
            println("--dcs: ${params.dcs}");
            throw new IllegalArgumentException('DCS/RCS sequence file must be specified')
        }
    }
    else if (params.library=='single_end')
    {
        params.pattern = '*_R1_*.fastq.gz'
    }
    else
    {
        println("--library: ${params.library}");
        throw new IllegalArgumentException('Library type must be either "paired_end", "single_end", "nanopore", or "single_cell"')
    }
    if (params.samples=='')
    {
        println("--samples: ${params.samples}");
        throw new IllegalArgumentException('Input sample directory must be specified')
    }
    if (params.metadata=='')
    {
        println("--metadata: ${params.metadata}");
        throw new IllegalArgumentException('Sample metadata file must be specified')
    }
    if (params.ref=='')
    {
        println("--ref: ${params.ref}");
        throw new IllegalArgumentException('CREATE reference directory must be specified')
    }
}
else if (params.discover)
{
    params.ref = ''
    if (params.ref=='')
    {
        println("--ref: ${params.ref}");
        throw new IllegalArgumentException('CREATE reference directory must be specified')
    }
    params.dcs = ''
    if (params.dcs=='')
    {
        println("--dcs: ${params.dcs}");
        throw new IllegalArgumentException('DCS/RCS sequence file must be specified')
    }
    params.long_reads = ''
    params.lr_pattern = '*.fastq.gz'
    if (params.long_reads=='')
    {
        println("--long_reads: ${params.long_reads}");
        throw new IllegalArgumentException('Long reads must be specified')
    }
    params.paired_reads = ''
    params.pr_pattern = '*_R{1,2}_*.fastq.gz'
    if (params.paired_reads=='')
    {
        println("--paired_reads: ${params.paired_reads}");
        throw new IllegalArgumentException('Paired reads must be specified')
    }
    params.prefixes = ''
    params.dump = ''
}
else
{
    throw new IllegalArgumentException('CREATE-seq mode must be specified either "--quant" or "--reference"')
}
if (params.reference && params.quant)
{
    println("--quant: ${params.quant}");
    println("--reference: ${params.reference}");
    throw new IllegalArgumentException('CREATE-seq mode must be specified either "--quant" XOR "--reference"')
}
if (!['OFF','INFO','DEBUG'].contains(params.log))
{
    println("--log: ${params.log}");
    throw new IllegalArgumentException('Log level must be one of: OFF, INFO, DEBUG')
}
if (params.outdir=='')
{
    println("--outdir: ${params.outdir}");
    throw new IllegalArgumentException('Output directory must be specified')
}


if (params.exec=="local")
{
    if (params.threads<0)
    {
        println("--threads: ${params.threads}");
        throw new IllegalArgumentException('Maximum number of threads must be at least 0')
    }
    if (params.memory<0)
    {
        println("--memory: ${params.memory}");
        throw new IllegalArgumentException('Memory must be at least 0 GB')
    }
    if (params.scratch!=false)
    {
        throw new IllegalArgumentException('Cannot assign scratch directory in local mode')
    }
    params.account = false
    if (params.account!='')
    {
        throw new IllegalArgumentException('Cannot assign account in local mode')
    }
    params.njobs = false
    if (params.njobs!=-1)
    {
        throw new IllegalArgumentException('Cannot assign number of jobs in local mode')
    }
}
else if (params.exec=="slurm")
{
    if (params.threads!=-1)
    {
        throw new IllegalArgumentException('Cannot assign thread limits in slurm mode')
    }
    if (params.memory!=-1)
    {
        throw new IllegalArgumentException('Cannot assign memory limits in slurm mode')
    }
}
else
{
    println("--exec: ${params.exec}");
    throw new IllegalArgumentException('Execution mode must be either "local" or "slurm"')
}
if (!['docker','conda'].contains(params.container))
{
    throw new IllegalArgumentException('Container must be either "docker" or "conda"')
}


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
params.each{k, v -> ['',false,-1,0].contains(v)?:print_val(k,v)}
log.info(" ")
params.manage_resources = (params.limits || params.exec!='local')


/*
 * Make CREATE workflow selector
 */
include { REFERENCE } from './workflow/reference/reference.nf'
include { QUANT } from './workflow/quant/quant.nf'
include { DISCOVER } from './workflow/discover/discover.nf'

workflow CREATE
{
    if (params.get('quant')!=null && params.get('quant')==true) QUANT()
    else if (params.get('reference')!=null && params.get('reference')==true) REFERENCE()
    else if (params.get('discover')!=null && params.get('discover')==true) DISCOVER()
}


/*
 * Run CREATE
 */
workflow { CREATE() }
