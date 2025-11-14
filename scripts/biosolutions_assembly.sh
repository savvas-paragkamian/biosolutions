#!/bin/bash -l

# Name: biosolutions_assembly.sh
# Purpuse: perform bacterial assembly to multiple isolates
# Author: Savvas Paragkamian
# Date: 01/04/2025
# Usage: ./biosolutions_assembly.sh assemblies_directories.txt 

######################### Data Validation ######################
#Download
#cd /media/sarlab/DATA
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

#dirs=$(find . -maxdepth 1 -type d ! -name "." | sed 's|^\./||')

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

# obtain the absolute path of the user supplied directories
dirs=$(realpath "$1")

# Check if the file exists
if [ ! -f "$dirs" ]; then
    echo "Error: File '$dirs' not found!"
    exit 1
fi

# Read the file line by line
#while IFS= read -r line; do
#    echo "$line"
#    echo $dirs
#done < "$dirs"

#exit 0
######################### Automated Short read Assembly ######################
cd /media/sarlab/DATA/biosolutions_project/genomes

# initiate conda in the script
source /opt/miniconda3/etc/profile.d/conda.sh

while IFS= read -r dir; do

	echo "Processing directory: $dir"

	conda activate perfect_assembly

	conda info

#	cd $dir
#	mkdir -p reads_qc
	# fastp quality
	echo "fastp of $dir"
#	fastp --in1 $dir/1.Cleandata/$dir.IS350_Clean.1.fq.gz \
#		--in2 $dir/1.Cleandata/$dir.IS350_Clean.2.fq.gz \
#		--out1 reads_qc/$dir.QC_1.fq.gz \
#		--out2 reads_qc/$dir.QC_2.fq.gz \
#		--unpaired1 reads_qc/$dir.QC_1_u.fq.gz \
#		--unpaired2 reads_qc/$dir.QC_2_u.fq.gz
#
#	# Unicycler for short-reads only:
#	unicycler -1 reads_qc/$dir.QC_1.fq.gz \
#		-2 reads_qc/$dir.QC_2.fq.gz \
#		-o unicycler_assembly
#
	conda activate quast
	conda info
#
#	#python /opt/miniconda3/envs/quast/bin/quast unicycler_assembly/assembly.fasta -t 12
#
#	cd ../
	
	echo "Finish Processing directory: $dir"
done < $dirs


echo "Finished All Assemblies"


# for the protein prediction use prodical
#source activate gtdbtk-2.3.2

#prodigal -i unicycler_assembly/assembly.fasta -o my.genes -a my.proteins.faa

######################### Taxonomy assignment ######################
# for the ani calculation

#gtdbtk ani_rep --genome_dir unicycler_assembly --out_dir ani_rep/ --cpus 15 -x fasta

# the file ani_rep/gtdbtk.ani_closest.tsv contains the taxonomy of the assembly

time_end=`date +%s`
time_exec=`expr $(( $time_end - $time_start ))`

echo "Execution time was $time_exec seconds"
