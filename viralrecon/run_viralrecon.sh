# This is the shell script used to submit jobs to UCL Myriad cluster for running the viralrecon pipeline on LUNAR & Usual Care samples  

#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=50G
#$ -pe smp 12
#$ -l tmpfs=50G
#$ -N viralrecon
#$ -wd /path/to/your/working/directory/Scratch/viralrecon_run       # Used the Scratch directory for this analysis

module load apptainer
source /path/to/your/miniforge3/etc/profile.d/conda.sh
conda activate nextflow

# Persistent caches so images and work survive between jobs
export NXF_SINGULARITY_CACHEDIR=/path/to/your/Scratch/singularity-cache
export APPTAINER_CACHEDIR=/path/to/your/Scratch/apptainer-cache
export APPTAINER_TMPDIR=$TMPDIR

# Launched from persistent Scratch, NOT $TMPDIR, so work/ and .nextflow/ survive
cd /path/to/your/working/directory/Scratch/viralrecon_run

nextflow run nf-core/viralrecon -r 3.0.0 -with-tower \
-profile singularity \
-params-file /path/to/your/data/directory/LUNAR_params.yml \
-c /path/to/your/data/directory/LUNAR.config \
--input /path/to/your/data/directory/LUNAR_PAN_samplesheet.csv \
--outdir /path/to/your/data/directory/NFresults \
-work-dir /path/to/your/working/directory/Scratch/viralrecon_run/work \
-resume
