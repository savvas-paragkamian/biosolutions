#!/bin/bash -l

# Name: biosolutions_assembly.sh
# Purpuse: perform bacterial assembly to multiple isolates
# Author: Savvas Paragkamian
# Date: 01/04/2025

# Usage: ./scripts/biosolutions_assembly.sh assemblies/sequences_directories.txt /media/sarlab/DATA/biosolutions_project/genomes/

######################### Data Validation ######################
#Download
#cd /to/save/the/DATA
#aws s3 cp  --recursive s3://bacmrjzd-598731762349/F25A910000023_BACmrjzD/ . 
#
#Validate download
#md5sum -c BGI_result/md5.txt 
#
#Uncompress directories
#cd BGI_result/Separate
#for file in *.tar.gz; do mkdir -p "${file%.tar.gz}" && tar -xvzf "$file" -C "${file%.tar.gz}"; done
#
#Check the number of sequencing files
#find . -type f -name "*.fq.gz"

######################### User Input ######################
time_start=`date +%s`

#cd /media/sarlab/DATA/
#md5sum -c BGI_result/md5.txt 
# if all is ok proceed to uncompress

# sparagkamian@sarrislab:/media/sarlab/DATA/biosolutions_project/genomes$ ls -d -1 */ > ~/Documents/biosolutions/assemblies/sequences_directories.txt 
# other way to create a list of directories
#dirs=$(find . -maxdepth 1 -type d ! -name "." | sed 's|^\./||')

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <filename> <path>"
    exit 1
fi

# Absolute paths
dirs=$(realpath "$1")
path=$(realpath "$2")

# Check if file exists
if [ ! -f "$dirs" ]; then
    echo "Error: File '$dirs' not found!"
    exit 1
fi

# Check if directory exists
if [ ! -d "$path" ]; then
    echo "Error: Directory '$path' not found!"
    exit 1
fi

# Prompt...

echo "User input is:"
echo "File of Directories of SRL microbes: $dirs"
echo "Path: $path"
microbes=$(wc -l < "$dirs")
echo "Number of SRL microbes: $microbes"

# Read the file line by line
#while IFS= read -r line; do
#    echo "$line"
#    #echo $dirs
#done < "$file"

######################### Automated Unicycler Assembly ######################
cd $path

# initiate conda in the script
source /opt/miniconda3/etc/profile.d/conda.sh
	
# activate the environment of assembly
conda activate autocycler
#exit 0

#----------------------------------------------------------------------#
#----------------------- Quality Filtering ----------------------------#
#----------------------------------------------------------------------#
while IFS= read -r dir_file; do

	dir=$(printf "%s" "$dir_file" | sed 's:[/[:space:]]*$::')
	echo "Processing directory: $dir"

	mkdir -p $dir/reads_qc

	# fastp quality filtering of short reads
	echo "fastp of $dir"
	fastp --in1 $dir/1.Cleandata/$dir*1.fq.gz \
		--in2 $dir/1.Cleandata/$dir*2.fq.gz \
		--out1 $dir/reads_qc/$dir.QC_1.fq.gz \
		--out2 $dir/reads_qc/$dir.QC_2.fq.gz \
		--unpaired1 $dir/reads_qc/$dir.QC_1_u.fq.gz \
		--unpaired2 $dir/reads_qc/$dir.QC_2_u.fq.gz \
		--thread 16

	# fastplong for filtering the long reads
	echo "fastplong of $dir"
	fastplong \
		-i $dir/1.Cleandata/$dir.filtered_reads.fq.gz \
		-o $dir/reads_qc/$dir.QC_long.fq.gz \
		--length_required 1000 \
		--qualified_quality_phred 10 \
		--thread 16

	echo "Finish Processing directory: $dir"

done < $dirs

#----------------------------------------------------------------------#
#----------------------- Unicycler Assembly ---------------------------#
#----------------------------------------------------------------------#
#cd $path

mkdir -p ../logs

rm -f ../logs/assemblies_unicycler_jobs.txt

# build the txt file with all the jobs
while IFS= read -r dir_file; do

	dir=$(printf "%s" "$dir_file" | sed 's:[/[:space:]]*$::')
	# Unicycler for short-reads only:
	echo "unicycler -1 $dir/reads_qc/$dir.QC_1.fq.gz -2 $dir/reads_qc/$dir.QC_2.fq.gz -l $dir/reads_qc/$dir.QC_long.fq.gz -o $dir/unicycler_assembly -t 8" >> ../logs/assemblies_unicycler_jobs.txt

done < $dirs

# run 3 jobs with GNU parallel
jobs=3

set +e
nice -n 19 parallel --jobs "$jobs" \
	--joblog ../logs/joblog.tsv \
	--results ../logs/logs \
	--timeout "$max_time" < ../logs/assemblies_unicycler_jobs.txt
set -e

#
	#conda activate quast
#
#	#python /opt/miniconda3/envs/quast/bin/quast unicycler_assembly/$dir_unicycler_assembly.fasta -t 12
#
#	cd ../
	
echo "Finished All Unicycler Assemblies"


time_end=`date +%s`
time_exec=`expr $(( $time_end - $time_start ))`

echo "Execution time was $time_exec seconds"
