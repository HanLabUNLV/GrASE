#' Create a Gencode TxDb object
#'
#' This function builds a transcript database (`TxDb`) object which you can then
#' use to build a Gencode `GenomicState` object. This function will download
#' the data from Gencode, import it into R, process it and build the `TxDb`
#' object.
#'
#' @param version A `character(1)` with the Gencode version number.
#' @param genome A `character(1)` with the human genome version number. Valid
#' options are `'hg38'` or `'hg19'`.
#' @param chrs A `character()` vector with the chromosome (contig) names to
#' keep.
#'
#' @return A [GenomicFeatures::TxDb-class] object.
#' @export
#' @author Leonardo Collado-Torres
#' @references Based on code for the `brainflowprobes` package at:
#' <https://github.com/LieberInstitute/brainflowprobes/blob/master/data-raw/create_sysdata.R>
#'
#' @examples
#' \dontrun{
#' ## Build from a local GTF file:
#' txdb <- gencode_txdb("path/to/gencode.v31.annotation.gtf.gz",
#'     genome = "hg19", chrs = "chr21")
#' txdb
#' }

gencode_txdb <- function( gtf_file, genome = c('hg19','hg38') , 
    chrs = paste0("chr", c(seq_len(22), "X", "Y", "M"))) {
    genome <- match.arg(genome)


    ## Import the data
    message(paste(Sys.time(), "importing", gtf_file))
    gencode_gtf <- rtracklayer::import(gtf_file)

    ## Keep only the main chrs
    message(paste(Sys.time(), "keeping relevant chromosomes"))
    gencode_gtf <- GenomeInfoDb::keepSeqlevels(gencode_gtf, chrs,
        pruning.mode = "coarse"
    )

    # Doesn't work because of the different seqlevels
    # txdb <- makeTxDbFromGFF(
    #     gtf_file,
    #     organism = 'Homo sapiens',
    #     chrominfo = Seqinfo(genome="hg19")
    # )

    message(paste(Sys.time(), "building the txdb object"))
    txdb <- txdbmaker::makeTxDbFromGRanges(gencode_gtf)
    return(txdb)
}


#' @export
#' @rdname gencode_txdb
#' @return A `character(1)` with the URL for the GTF Gencode file of interest.
#'
#' @examples
#'
#' ## Locate the GTF file for Gencode version 31 for hg19
#' gencode_source_url(version = "31", genome = "hg19")
gencode_source_url <- function(version = "34", genome = c("hg38", "hg19")) {
    genome <- match.arg(genome)
    source_url <- paste0(
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
        "release_", version, if (genome == "hg19") "/GRCh37_mapping",
        "/gencode.v", version, if (genome == "hg19") "lift37",
        ".annotation.gtf.gz"
    )

    return(source_url)
}


