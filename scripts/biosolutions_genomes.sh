#!/bin/bash -l

# Name: biosolutions_genomes.sh
# Purpuse: perform genome statistics, quality and annotation 
# with different tools
# Author: Savvas Paragkamian
# Date: 23/11/2025

# Usage: ./scripts/biosolutions_genomes.sh biosolutions_batchfile.txt

######################### User Input ######################
time_start=`date +%s`

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <batchfile/with/assembly/paths>"
    exit 1
fi

# Absolute paths
assemblies=$(realpath "$1")

# Check if file exists
if [ ! -f "$assemblies" ]; then
    echo "Error: File '$assemblies' not found!"
    exit 1
fi

# Prompt...

echo "User input is:"
echo "Batchfile of SRL microbes: $assemblies"
microbes=$(wc -l < "$assemblies")
echo "Number of SRL microbes: $microbes"


# initiate conda in the script
source /opt/miniconda3/etc/profile.d/conda.sh
	
# activate the environment of assembly

# --------------------- Assembly statistics --------------------#
# quast takes about a second per assembly

#conda activate quast
#
#while IFS= read -r dir_file; do
#
#	echo $dir_file
#
#	dir=$(echo $dir_file | awk '{print $1}')
#
#	echo "Processing assembly: $dir"
#	dir_path=$(echo $dir | awk -F"/" '{print $1 "/" $2}')
#	echo "Directory of assembly: $dir_path"
#
#	cd $dir_path
#
#	assembly=$(printf '%s' "$dir" | sed "s|$dir_path/||g")
#
#	echo "assembly: $assembly"
#
#	mkdir -p quast
#
#	# quast
#	echo "quast analysis of $dir"
#	python /opt/miniconda3/envs/quast/bin/quast \
#		$assembly \
#		-t 12 \
#		-o quast 
#
#	cd ../../
#
#	echo "Finish Processing directory: $dir_path"
#
#done < $assemblies
#
#
#	cd ../

# ------------------------ Annotation ----------------------------#
########################## With Bakta #############################

#conda activate bakta
#
#while IFS= read -r dir_file; do
#
#	echo $dir_file
#
#	dir=$(echo $dir_file | awk '{print $1}')
#
#	echo "Processing assembly: $dir"
#	dir_path=$(echo $dir | awk -F"/" '{print $1 "/" $2}')
#	echo "Directory of assembly: $dir_path"
#
#	cd $dir_path
#
#	assembly=$(printf '%s' "$dir" | sed "s|$dir_path/||g")
#
#	echo "assembly: $assembly"
#
#	# bakta
#	echo "bakta analysis of $dir"
#	bakta $assembly  \
#		--db /media/sarlab/DATA/databases/bakta_v6.0/db \
#		--threads 12 \
#		--output bakta \
#		--force
#
#	cd ../../
#
#	echo "Finish Processing directory: $dir_path"
#
#done < $assemblies

########### merge annotations to one file #####################

out="bakta_merged.tsv"
first=1
> "$out"

while IFS= read -r dir_file; do

	echo $dir_file

	dir=$(echo $dir_file | awk '{print $1}')

	echo "Processing assembly: $dir"
	dir_path=$(echo $dir | awk -F"/" '{print $1 "/" $2}')
	echo "Directory of assembly: $dir_path"

	annotationofassembly="$dir_path/bakta/assembly.tsv"
	assembly=$(printf '%s' "$dir_path" | sed "s|"genomes/"||g")

	echo "annotations of assembly : $annotationofassembly"

	# bakta
	echo "bakta annotation of $assembly"
  	# safety: skip if missing/empty
  	if [ ! -s "$annotationofassembly" ]; then
    		echo "WARNING: missing/empty: $annotationofassembly" >&2
    		continue
  	fi

	# handle differently the first for the header
	if [ "$first" -eq 1 ]; then
    		# keep the 6-line header ONCE, but prepend "assembly" to the
		#header line that holds column names
    		# the 6th line is the TSV column header line.
    		awk -v a="$assembly" -v OFS='\t' '
		NR==6 { print "assembly", $0; next }
		NR<=5 { print; next }
		{ next }
		' "$annotationofassembly" >> "$out"
		first=0
	fi
	

	# append data rows (skip first 6 lines), prepend assembly column
	tail -n +7 "$annotationofassembly" \
    		| awk -v a="$assembly" -v OFS='\t' '{print a, $0}' >> "$out"

  	echo "Finish Processing directory: $dir_path"

done < $assemblies


########################## CheckM2 #############################

#conda activate checkm2
#
#while IFS= read -r dir_file; do
#
#	echo $dir_file
#
#	dir=$(echo $dir_file | awk '{print $1}')
#
#	echo "Processing assembly: $dir"
#	dir_path=$(echo $dir | awk -F"/" '{print $1 "/" $2}')
#	echo "Directory of assembly: $dir_path"
#
#	cd $dir_path
#
#	assembly=$(printf '%s' "$dir" | sed "s|$dir_path/||g")
#
#	echo "assembly: $assembly"
#
#	# bakta
#	echo "checkM2 analysis of $dir"
#	checkm2 predict \
#		--threads 12 \
#		--input $assembly \
#		--output-directory checkm2 \
#		--force
#
#	cd ../../
#
#	echo "Finish Processing directory: $dir_path"
#
#done < $assemblies
#
########################## antiSMASH #############################

#conda activate antismash
#
#while IFS= read -r dir_file; do
#
#	echo $dir_file
#
#	dir=$(echo $dir_file | awk '{print $1}')
#
#	echo "Processing assembly: $dir"
#	dir_path=$(echo $dir | awk -F"/" '{print $1 "/" $2}')
#	echo "Directory of assembly: $dir_path"
#
#	cd $dir_path
#
#	#assembly=$(printf '%s' "$dir" | sed "s|$dir_path/||g")
#
#	#echo "assembly: $assembly"
#
#	# bakta
#	echo "antismash analysis of $dir"
#	antismash bakta/assembly.embl  \
#		--cpus 12 \
#		--output-dir antismash \
#		--cc-mibig \
#		--cb-general \
#		--enable-genefunctions \
#		--cb-subclusters \
#		--cb-knownclusters \
#		--enable-lanthipeptides \
#                --enable-lassopeptides \
#		--enable-nrps-pks \
#		--enable-sactipeptides \
#		--enable-t2pks \
#                --enable-thiopeptides \
#		--enable-tta 
##		--fullhmmer \
#
#	cd ../../
#
#	echo "Finish Processing directory: $dir_path"
#
#done < $assemblies

########################## ending #############################
time_end=`date +%s`
time_exec=`expr $(( $time_end - $time_start ))`

echo "Execution time was $time_exec seconds"
