<!--
  REQUIRED NOTICE: Copyright (c) 2020-2023, Regents of the University of California
  All rights reserved. https://polyformproject.org/licenses/noncommercial/1.0.0
  
  This software was developed by the Daniel Kim lab at the University of California, Santa Cruz.
  Authors: Roman E. Reggiardo, Vikas Peddu, Alex D. Hill
  
  The licensor grants you a copyright license for the software to do everything you might do with
  the software that would otherwise infringe the licensor’s copyright in it for any permitted
  purpose.
  
  As far as the law allows, the software comes as is, without any warranty or condition, and the
  licensor will not be liable to you for any damages arising out of these terms or the use or
  nature of the software, under any kind of legal claim.
-->


# `CREATE`-seq
## Comprehensive REptitive And Transposable Element Sequencing
### Alex D. Hill

CREATE is a pipeline for the quantification of repetitive elements alongside canonical genes. It is designed to be very high throughput, and can be run very quickly on a large number of samples. The core of CREATE's functions are based on Trim-Galore, Salmon, and Tximeta for fast read processing and storage.

There are two modes to CREATE:
  - REFERENCE: This mode creates a rich 'complete' set of reference files that contain a formatted compilation of both GENCODE annotations and transcripts, as well as repetitive element annotations and transcripts pulled from the RepeatMasker track on the UCSC genome browser. Reference directories can also be downloaded for ease of use.
  - QUANT: This mode takes a set of raw reads and outputs a highly efficient H5+RDS file that efficiently store a SummarizedExperimet (or SingleCellExperiment) object that can be used for downstream analysis.

#### Installation:
As a nextflow pipeline, CREATE can be run with only an installation of conda. Make an environment with Nextflow installed and force create to utilize conda and the container format. *NOTE*: a docker installation is highly recommended to guarantee compatibility with most systems.

#### examples: REFERENCE
```
nextflow run alexdhill/CREATE --reference \
  --genome HG38 \
  --version 39 \
  --index short \
  --outdir HG38_reference.crt
```

#### examples: QUANT
```
nextflow run alexdhill/CREATE --quant \
  -resume \
  --samples <path/to/samples> \
  --ref <path/to/reference> \
  --library paired_end \
  --outdir create_results
```


#### usage:
<!--    discover          Discover novel isoforms [BETA] -->
```
 $ nextflow run alexdhill/CREATE --help
Usage: nextflow run alexdhill/CREATE [--quant/--reference/--analyze] [options]
    reference         Build a reference for CREATE
    quant             Quantify RNA-seq reads with CREATE

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

    --quant [default]
        --samples     Input samples directory [REQUIRED]
        --ref         CREATE reference directory [REQUIRED]
        --library     Library type [paired_end,nanopore,single_cell] [REQUIRED]
                      NOTE: the reference provided must contain a compatible index
        --pattern     File pattern for input samples
        --metadata    An unnamed, two-column CSV file of samples,condition

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
```

#### roadmap:
- [X] Finish Reference workflows
  - [X] Finish salmon reference for short reads
  - [X] Finish minimap reference for long reads
  - [X] Finish splintr reference for single-cell
- [ ] Finish Quant workflows
  - [X] Finish Paired read quant
  - [X] Finish Single read quant
  - [X] Finish Single-Cell quant
  - [X] Finish Nanopore quant
    - [X] Add DCS/RCS removal
- [ ] Add novel isoform discovery (FLAIR)
- [X] Add Analysis ~~workflow~~ step to quant
  - [X] Finish biotype plots
  - [ ] Finish PCA plots
  - [X] Finish Volcano plots
- [X] Flesh README
- [ ] ~~Add profile selection to nextflow config~~ Add read-based resource requirements to processes
