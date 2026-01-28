#!/bin/bash
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G   # memory pool for all cores
#SBATCH --time=7-00:00:00   # time limit
#SBATCH --mail-user=clarayu@caltech.edu   # email address
#SBATCH --mail-type=FAIL   # Notify on failure.
#SBATCH --output=/resnick/groups/MazmanianLab/clarayu/slurm-%j.out

export PATH="/resnick/groups/MazmanianLab/clarayu/miniforge3/bin:$PATH"
source /home/${USER}/.bashrc
source /resnick/groups/MazmanianLab/clarayu/miniforge3/etc/profile.d/conda.sh
conda activate snakemake
which conda
conda info --json > /dev/null && echo "Conda check passed" || echo "Conda check failed"
export SNAKEMAKE_OUTPUT_CACHE=/resnick/groups/MazmanianLab/clarayu/microbiome_sample_search/.snakemake_cache

cd /resnick/groups/MazmanianLab/clarayu/microbiome_sample_search/metagen-paladin
snakemake --unlock
snakemake \
  --use-conda \
  --conda-frontend conda \
  --conda-prefix $HOME/.snakemake-conda \
  --jobs 9999 --latency-wait 300 --rerun-incomplete --resources load=1000 --keep-going \
  --cluster-config config/cluster.json \
  --cluster "sbatch -p {cluster.queue} --time={cluster.time} --nodes=1 --cpus-per-task={cluster.cores} --mem={cluster.mem} -J {cluster.jobname} -o {cluster.outfile} -e {cluster.errfile}" \
  --printshellcmds --reason --verbose
  #--cluster "sbatch -p {cluster.queue} --time={cluster.time} --mem={cluster.mem} --ntasks={cluster.cores} -J {cluster.jobname} -o {cluster.outfile} -e {cluster.errfile}"

#snakemake --use-conda -k --cluster "sbatch -p expansion --time=6:00:00 --mem=100GB --cpus-per-task=4" --jobs 9999 --rerun-incomplete

# snakemake --use-conda -k --cluster "sbatch -p expansion --time=8:00:00 --mem=200GB --cpus-per-task=8" --jobs 9999 --rerun-incomplete
