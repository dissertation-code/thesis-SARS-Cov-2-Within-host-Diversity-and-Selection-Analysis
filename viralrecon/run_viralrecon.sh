#!/bin/bash -l

#$ -l h_rt=12:00:00
#$ -l mem=100G
#$ -pe smp 12
#$ -l tmpfs=100G
#$ -N viralrecon
#$ -wd /path/to/your/working/directory/ # Used the Scratch directory

cd $TMPDIR

module load apptainer

source /path/to/your/miniforge3/etc/profile.d/conda.sh

conda activate nextflow

export APPTAINER_CACHEDIR=/path/to/your/apptainer-cache  # Used the Scratch directory 
export APPTAINER_TMPDIR=/path/to/your/apptainer-tmp      # Used the Scratch directory 

nextflow run nf-core/viralrecon -r 3.0.0 -with-tower \
-profile singularity \
-params-file /path/to/your/data/directory/LUNAR_params.yml \
-c /path/to/your/data/directory/LUNAR.config \
--input /path/to/your/data/directory/LUNAR_PAN_samplesheet.csv \
--outdir /path/to/your/data/directory/NFresults
