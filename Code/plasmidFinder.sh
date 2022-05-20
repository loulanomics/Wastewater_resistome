#!/usr/bin/env bash

##############################
### align contigs to 
### PlasmidFinder database
### Lou LaMartina, May 19 2022
##############################


db=Data/PlasmidFinder/plasmid_database.fsa


# make blast database with plasmid sequences
makeblastdb -dbtype nucl -in $db -input_type fasta


# # # # # # # 


in=Data/Contigs
out=Data/output_plasmidFinder


for q in $in/*.fa; do
	f="$(basename $q)"
	echo "Aligning to $f . . ."
	blastn -db $db -query $q -task blastn -perc_identity 90 -outfmt 6 -out $out/plasmids_${f%.*}.txt
	echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n$(cat $out/plasmids_${f%.*}.txt)" > $out/plasmids_${f%.*}.txt
done


# remove blastn artifacts
rm Data/PlasmidFinder/*.fasta.n*

