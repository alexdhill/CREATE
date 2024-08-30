main = function()
{
    source(args$helpers)
    suppressMessages({
        library(argparse)
    })

    parser = ArgumentParser()
    parser$add_argument("--quants", help = "The `CREATE quant` dir", type = "character")
    parser$add_argument("--reference", help = "The CREATE reference dir", type = "character")
    parser$add_argument("--helpers", help = "The Kim Lab R helpers file", type = "character")
    parser$add_argument("--output", help = "The output directory", type = "character")
    args = parser$parse_args()

    message("Importing quantifications...")
    quants = loadQuantH5(args$quants)

    message("Importing transcript to gene and biotype map")
    biotypes = readTranscriptMap(args$reference)

    message("Running differential expression with DESeq2...")
    
}

main()