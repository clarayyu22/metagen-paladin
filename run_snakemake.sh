#!/bin/bash
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G   # memory pool for all cores
#SBATCH --time=7-00:00:00   # time limit
#SBATCH --mail-user=clarayu@caltech.edu   # email address
#SBATCH --mail-type=FAIL   # Notify on failure.
#SBATCH --output=/dev/null

source /home/${USER}/.bashrc
source activate snakemake

export SNAKEMAKE_OUTPUT_CACHE=/central/groups/MazmanianLab/clara/microbiome_sample_search/.snakemake_cache

cd /central/groups/MazmanianLab/clara/microbiome_sample_search/metagen-paladin
snakemake --unlock
snakemake --use-conda -k --jobs 9999 --latency-wait 300 --rerun-incomplete --resources load=100 \
  --cluster-config config/cluster.json \
  --cluster "sbatch -p {cluster.queue} --time={cluster.time} --nodes=1 --cpus-per-task={cluster.cores} --mem={cluster.mem} -J {cluster.jobname} -o {cluster.outfile} -e {cluster.errfile}"
  --printshellcmds --reason --verbose
  #--cluster "sbatch -p {cluster.queue} --time={cluster.time} --mem={cluster.mem} --ntasks={cluster.cores} -J {cluster.jobname} -o {cluster.outfile} -e {cluster.errfile}"

#snakemake --use-conda -k --cluster "sbatch -p expansion --time=6:00:00 --mem=100GB --cpus-per-task=4" --jobs 9999 --rerun-incomplete

# snakemake --use-conda -k --cluster "sbatch -p expansion --time=8:00:00 --mem=200GB --cpus-per-task=8" --jobs 9999 --rerun-incomplete