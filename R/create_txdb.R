## create_txdb.R -- create txdb SQLite file
##
## These functions are private; this is just for reproducibility and updating
## this package.

library(GenomicFeatures)
library(rtracklayer)

.makeTxDb <- function(gtf, len_file, data_source, species) {
	chrominfo <- read.delim(len_file, header=FALSE, col.names=c("chrom", "length", "is_circular"),
													colClasses=c("character", "integer", "logical"))

	# this code adapted from makeTranscriptDbFromGFF, but adapted to subset by protein coding
	exonRankAttributeName <- "exon_number"

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

# load GTF
# These features only exclude start and stop codons in the GTF file.
feature.type <- c("gene", "mRNA", "exon", "CDS")
gtf <- import("inst/extdata/Zea_mays.AGPv2.17.gtf.gz", format="gtf", feature.type=feature.type, asRangedData=FALSE)

source_url <- "ftp://ftp.ensemblgenomes.org/pub/plants/release-17/gtf/zea_mays/Zea_mays.AGPv2.17.gtf.gz"
# create txdb from protein coding sequences from ensembl
message("creating ensemble whole gene set txdb...")
txdb <- .makeTxDb(gtf, "inst/extdata/Zea_mays.AGPv2.17.lengths.txt",
									source_url,
									"Zea mays")
message("saving ensemble whole gene set txdb...")
saveDb(txdb, "inst/extdata/TxDb.Zmays.Ensembl.AGPv2.17.sqlite")

# using the filtered gene set ID's from maizesequence.org
# (http://ftp.maizesequence.org/release-5b/filtered-set/ZmB73_5b_FGS.gff.gz)
# create another txdb object of filtered genes
# the following was used to create the filtered gene ID list
# gzcat ZmB73_5b_FGS.gff.gz |  grep gene | cut -f9 | cut -f1 -d";" | sed 's/ID=//' | sort | uniq > ZmB73_5b_FGS_gene_ids.txt
message("creating ensemble filtered gene set txdb using filtered gene IDs from maizesequence.org...")

filtered_gene_ids <- scan("inst/extdata/ZmB73_5b_FGS_gene_ids.txt", what=character())

txdb_filtered <- .makeTxDb(gtf[gtf$gene_id %in% filtered_gene_ids],
									"inst/extdata/Zea_mays.AGPv2.17.lengths.txt",
									source_url,
									"Zea mays")
message("saving ensemble filtered gene set txdb...")
saveDb(txdb_filtered, "inst/extdata/TxDb.Zmays.EnsemblFiltered.AGPv2.17.sqlite")



