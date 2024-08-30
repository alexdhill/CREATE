suppressPackageStartupMessages({
    library(tidyr)
    library(dplyr)
    library(stringr)
    library(argparse)
    library(txdbmaker)
    library(GenomicFeatures)
    library(BiocFileCache)
})

main <- function()
{
    parser = ArgumentParser()
    parser$add_argument(
        "-a", "--annotation", action="store",
        help="The reference annotation"
    )
    args = parser$parse_args()

    build = args$annotation %>%
        strsplit("_") %>%
        unlist() %>%
        head(n=1)

    txdb = txdbmaker::makeTxDbFromGFF(
        args$annotation,
        dataSource="de-novo",
        organism="Homo sapiens"
    )

    dir.create(".bioccache")
    bfc = BiocFileCache(".bfccache")
    saveDb(txdb, file=bfcnew(bfc, rname=annotation, ext=".sqlite"))
}

main()