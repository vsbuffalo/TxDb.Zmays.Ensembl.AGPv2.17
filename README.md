# TxDb.Zmays.Ensembl.AGPv2.17

This TranscriptDb was generated with AGPv2.17 GTF from Ensembl, downloaded with:

    wget -O - ftp://ftp.ensemblgenomes.org/pub/plants/release-17/gtf/zea_mays/Zea_mays.AGPv2.17.gtf.gz \
      > inst/extdata/Zea_mays.AGPv2.17.gtf

This GTF was downloaded on 2014-03-08. This version is the last release from
Ensemble for AGPv2.

Chromosome length information (`inst/extdata/AGPv2.17.lengths.txt`) was made
with:

    wget -O - ftp://ftp.ensemblgenomes.org/pub/plants/release-17/fasta/zea_mays/dna/Zea_mays.AGPv2.17.dna.toplevel.fa.gz \
      | bioawk -c fastx '{print $name"\t"length($seq)"\tFALSE"}' > inst/extdata/Zea_mays.AGPv2.17.lengths.txt

## Filtered and Whole Gene Sets

[Maizesequence.org](www.maizesequence.org) (now
[Gramene](ensembl.gramene.org/Zea_mays/Info/Index)) has different annotation
than Ensembl, including classification of genes into filtered whole gene sets.
Currently, there are `63,331` protein-coding genes in Ensembl's data set (see
with `length(genes(txdb))`). Only a subset of these genes are supported by
evidence; others may be pseudogenes or gene fragments (see [this description on
CoGe](http://genomevolution.org/wiki/index.php/Sequenced_plant_genomes#Maize.2FCorn)
for more information). The file `inst/extdata/ZmB73_5b_FGS.gff.gz` from
MaizeSequence.org was processed so as to gather all filtered gene IDs:

    gzcat ZmB73_5b_FGS.gff.gz | grep gene | cut -f9 | \
      cut -f1 -d";" | sed 's/ID=//' | sort | uniq > ZmB73_5b_FGS_gene_ids.txt

The file `ZmB73_5b_FGS.gff.gz` was downloaded on 2014-03-18.

With these filtered IDs, a second TranscriptDb of filtered genes is created:
`TxDb.Zmays.EnsemblFiltered.AGPv2.17`.

## Using this Package

Install with:

    R CMD INSTALL TxDb.Zmays.Ensembl.AGPv2.17

Then, to use:

    library(TxDb.Zmays.Ensembl.AGPv2.17)
    txdb <- TxDb.Zmays.Ensembl.AGPv2.17
    # or
    # txdb <- TxDb.Zmays.EnsemblFiltered.AGPv2.17
    exons <- exonsBy(txdb, "gene")
