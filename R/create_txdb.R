## create_txdb.R -- create txdb SQLite file
##
## These functions are private; this is just for reproducibility and updating
## this package.

library(GenomicFeatures)
library(rtracklayer)

.makeTxDb <- function(gtf_file, len_file, data_source, species) {
	chrominfo <- read.delim(len_file, header=FALSE, col.names=c("chrom", "length", "is_circular"),
													colClasses=c("character", "integer", "logical"))

	# this code adapted from makeTranscriptDbFromGFF, but adapted to subset by protein coding
	exonRankAttributeName <- "exon_number"
	feature.type <- c("gene", "mRNA", "exon", "CDS")
  gtf <- import(gtf_file, format="gtf", feature.type=feature.type, asRangedData=FALSE)

	protein_coding <- gtf[gtf$source == "protein_coding"]
	if (all(c("gene_id", "transcript_id") %in% colnames(mcols(protein_coding)))) {
		tables <- GenomicFeatures:::.prepareGTFTables(protein_coding, exonRankAttributeName)
	}

	metadata <- GenomicFeatures:::.prepareGFFMetadata(gtf_file, data_source, species, NA)
	if (is.na(chrominfo)) {
		message("Now generating chrominfo from available sequence names. No chromosome length information is available.")
		chroms <- unique(tables[["transcripts"]][["tx_chrom"]])
		chrominfo <- data.frame(chrom=chroms, length=rep(NA,
																										 length(chroms)), is_circular=GenomicFeatures:::matchCircularity(chroms,
																										 GenomicFeatures:::DEFAULT_CIRC_SEQS))
	}
	txdb <- makeTranscriptDb(transcripts=tables[["transcripts"]],
													 splicings=tables[["splicings"]], genes=tables[["genes"]],
													 chrominfo=chrominfo, metadata=metadata, reassign.ids=TRUE)
	return(txdb)
}

message("creating txdb...")
txdb <- .makeTxDb("inst/extdata/Zea_mays.AGPv2.17.gtf", "inst/extdata/Zea_mays.AGPv2.17.lengths.txt",
									"ftp://ftp.ensemblgenomes.org/pub/plants/release-17/gtf/zea_mays/Zea_mays.AGPv2.17.gtf.gz",
									"Zea mays")

message("saving txdb...")
saveDb(txdb, "inst/extdata/TxDb.Zmays.Ensembl.AGPv2.17.sqlite")
