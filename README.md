# MetaGen-Fetch + Paladin - Processing public metagenomes


----
**Give credit to source**

** rename config_template.yaml to config.yaml and udpdate values before execution


This is a Snakemake workflow for downloading and quality-controlling metagenomic read files (FASTQ) from ENA. It uses [fastq-dl](https://github.com/rpetit3/fastq-dl) to first download a set of paired-end reads from ENA, and subsequently QCs the data using the [metagen-fastqc](https://github.com/alexmsalmeida/Metagen-FastQC) script.

## Installation

1. Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (tested v7.32.3)

2. Prepare and index a FASTA file of your host genome for decontamination. You can follow the instructions provided in the metagen-fastqc script [here](https://github.com/alexmsalmeida/Metagen-FastQC).

3. Clone this repository
```
git clone https://github.com/alexmsalmeida/metagen-fetch.git
```

## How to run

1. Edit the configuration file [`config/config.yml`](config/config.yml).
    - `input_file`: TSV file (no header) with run accessions listed as the first column and corresponding study accessions as the second column.
    - `output_dir`: Output directory to store cleaned files.
    - `host_ref`: Location of the indexed FASTA file for host decontamination.

2. Run the pipeline on a cluster (e.g., SLURM)
```
snakemake --use-conda -k -j 25 --profile config/slurm
```

3. Cleaned files will be stored in the specified output directory followed by `[study]/[run]`.
