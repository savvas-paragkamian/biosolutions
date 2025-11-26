#!/bin/bash

## Purpose of script: 
## Run the gtdb tk for the taxonomic inference and 
## de nove tree estimation.
##
## Author: Savvas Paragkamian

# ------------------ Batchfile -----------------#
# create a batch file for the assemblies 
# inside the directory of the project run the following command 
# to wirte to a batch file the path of the assembly and the id of 
# the microbe
#find -L . -iname "*assembly.fasta" | sort | grep "genomes" | gawk -F"/" '{gsub("\\./","",$0); print $0 "\t" $2}' > biosolutions_batchfile.txt
# the -L is to follow symbolic links also

# ------------------ de novo wf ----------------#
# this run takes 10 hours

gtdbtk de_novo_wf \
	--extension fasta \
	--batchfile biosolutions_batchfile.txt \
	--out_dir gtdbtk_denovo \
	--bacteria \
	--outgroup_taxon p__Chloroflexota \
	--cpus 16 

# outgroup is hard to define https://github.com/Ecogenomics/GTDBTk/issues/390

# ------------------ classify ----------------#
# this run takes 1 hour

#gtdbtk classify_wf \
#	--batchfile biosolutions_batchfile.txt \
#	--out_dir gtdbtk_classify \
#	--cpus 16

