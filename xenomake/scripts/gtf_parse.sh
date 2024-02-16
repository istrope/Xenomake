#!/bin/bash
#WARNING: This assumes that for "gene" features, the third attribute is the gene name (if present)
# If the third attribute is not the gene name, then the gene id (first column) is used.
#Extracts Gene Features from GTF file and returns in tabular format... used for converting between ENSEMBL and GENE SYMBOLS
cat $@ | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' |sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_type "//'| sed 's/gene_name "//' | sed 's/"//g' | awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' | sed "1i\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength"
