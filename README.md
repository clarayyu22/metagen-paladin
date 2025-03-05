# MetaGen-Fetch + Paladin - Processing and aligning public metagenomes


----
Credits: alexmsalmeida (metagen-fetch), ToniWestbrook (paladin)

This is a Snakemake workflow for downloading, quality-controlling metagenomic read files (FASTQ) from ENA, and using PALADIN to align them to FASTA files. It uses [fastq-dl](https://github.com/rpetit3/fastq-dl) to first download a set of paired-end reads from ENA, and subsequently QCs the data using the [metagen-fastqc](https://github.com/alexmsalmeida/Metagen-FastQC) script. It then uses [paladin] (https://github.com/ToniWestbrook/paladin) to align the cleaned metagenomic files to prepared biological sequences. 

## Installation

1. Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (tested v7.32.3)

2. Prepare and index a FASTA file of your host genome for decontamination. You can follow the instructions provided in the metagen-fastqc script [here](https://github.com/alexmsalmeida/Metagen-FastQC).

3. Clone this repository
```
git clone https://github.com/clarayyu22/metagen-paladin.git
```

4. Install PALADIN inside the metagen-paladin directory. Follow the instructions prvoided [here] (https://github.com/ToniWestbrook/paladin). 

## How to run

1. Edit the configuration file [`config_template/config.yml`](config_template/config.yml).
    - Rename `config_template/config.yml` to `config/config.yml`
    - `input_file`: TSV file (no header) with run accessions listed as the first column and corresponding study accessions as the second column.
    - `output_fastq_dir`: Output directory to store cleaned files.
    - `host_ref`: Location of the indexed FASTA file for host decontamination.
    - `fasta_dir`: Home folder of the directory with the FASTA file containing sequences of interest. 
    - `output_align_dir`: Output directory to store alignment files.
    - `fasta`: Name of FASTA file (text before .fasta). 
    - Directory paths should end in "`\`"

2. Prepare the FASTA file containing sequences of interest
    - Place the FASTA file inside a directory named `fasta` inside of `fasta_dir`
        - Path format: `fasta_dir` `fasta`/`fasta`.fasta.txt
    - Run the following command: 
        ```
        ./paladin/paladin prepare -r1 -f [FASTA file path]`
        ```

3. Run the pipeline on a cluster (e.g., SLURM)
    - Activate snakemake
    - Run the following command: 
        ```
        snakemake --use-conda -k -j 25 --profile config/slurm
        ```

3. Cleaned files will be stored in the specified FASTQ output directory followed by `[study]/[run]`. FASTQ files will be deleted after alignment to save memory. Alignment output files will be stored in the specified align output directory and named `[fasta]_[fastq sequence]_out_uniprot.tsv` and `[fasta]_[fastq sequence]_out.sam`. 
