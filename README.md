# TxDb.Zmays.Ensembl.AGPv2.17

This txdb was generated with AGPv2.17 GTF from Ensembl, downloaded with:

    wget -O - ftp://ftp.ensemblgenomes.org/pub/plants/release-17/gtf/zea_mays/Zea_mays.AGPv2.17.gtf.gz \
      > inst/extdata/Zea_mays.AGPv2.17.gtf

This GTF was downloaded on 2014-03-08. This version is the last release from
Ensemble for AGPv2.

Chromosome length information (`inst/extdata/AGPv2.17.lengths.txt`) was made
with:

    wget -O - ftp://ftp.ensemblgenomes.org/pub/plants/release-17/fasta/zea_mays/dna/Zea_mays.AGPv2.17.dna.toplevel.fa.gz \
      | bioawk -c fastx '{print $name"\t"length($seq)"\tFALSE"}' > inst/extdata/Zea_mays.AGPv2.17.lengths.txt
