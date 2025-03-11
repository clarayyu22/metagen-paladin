import os

configfile: 'config/config.yml'
configfile: 'config/cluster.json'

cpus_ena = config['ena_download']['cores']
cpus_qc = config['metagen_qc']['cores']
cpus_paladin = config['paladin']['cores']
host_ref = config['host_ref']

INPUT_FILE = config['input_file']
OUTPUT_FASTQ_DIR = config['output_fastq_dir']
PALADIN_REF_DIR = config['paladin_ref_dir']
FASTA_PATH = config['fasta_path']
FASTA_NAME = os.path.basename(FASTA_PATH)
config["paladin_db"] = os.path.join(config['paladin_ref_dir'], FASTA_NAME)

SAMPLES = []
STUDIES = []
samp2study = {}

os.system("chmod -R +x scripts")

with open(INPUT_FILE) as f:
    for line in f:
        cols = line.rstrip().split("\t")
        samp2study[cols[0]] = cols[1]
        SAMPLES.append(cols[0])
        STUDIES.append(cols[1])

for sample in samp2study:
    dirname = OUTPUT_FASTQ_DIR+samp2study[sample]+"/"+sample+"/logs"
    if not os.path.exists(dirname):
        os.makedirs(dirname)    

rule targets:
    input:
        expand([OUTPUT_FASTQ_DIR+"{study}/{sample}/paladin_done.txt"], zip, study=STUDIES, sample=SAMPLES)
   
rule ena_download:
    output:
        fwd = "{output}/{study}/{sample}/{sample}_1.fastq.gz"
    params:
        outdir = "{output}/{study}/{sample}",
        sample = "{sample}"
    conda:
        "config/envs/ena_download.yml"
    resources:
        ncores = cpus_ena
    shell:
        """
        # Download files
        fastq-dl --cpus {resources.ncores} -a {wildcards.sample} --outdir {params.outdir}
        
        # Check if SRA file exists and process it if found
        if [ -f "{params.outdir}/{params.sample}.sra" ]; then
            cd {params.outdir}
            fastq-dump --split-3 *.sra --gzip
            rm {params.outdir}/{params.sample}.sra
        fi
        
        # Define possible file patterns
        R1="{params.outdir}/{params.sample}_1.fastq.gz"
        R2="{params.outdir}/{params.sample}_2.fastq.gz"
        SINGLE="{params.outdir}/{params.sample}.fastq.gz"
        
        # Immediately remove R2 if it exists
        if [ -f "$R2" ]; then
            rm -f "$R2"
        fi
        
        # If both R1 and singlet exist, use R1
        if [ -f "$R1" ] && [ -f "$SINGLE" ]; then
            rm -f "$SINGLE"
            
        # If only singlet exists, use that
        elif [ -f "$SINGLE" ]; then
            cp "$SINGLE" {output.fwd}
        
        # Scenario 2: If only R1 exists, use R1
        elif [ -f "$R1" ]; then
            # No need to do anything, R1 is already named properly
            :

        else
            echo "No suitable forward read file found" >&2
            exit 1
        fi
        
        # Clean up original files after successful copy
        if [ -f "$SINGLE" ]; then
            rm -f "$SINGLE"
        fi
        """

rule metagen_qc:
    input:
        fwd = "{output}/{study}/{sample}/{sample}_1.fastq.gz"
    output:
        msg = "{output}/{study}/{sample}/qc_done.txt"
    params:
        bwa_ref = host_ref,
        fwd = "{output}/{study}/{sample}/{sample}_1_clean.fastq.gz"
    conda:
        "config/envs/metagen_qc.yml"
    resources:
        ncores = cpus_qc
    shell:
        """
        ./scripts/metagen-fastqc.sh -t {resources.ncores} -c {params.bwa_ref} -f {input.fwd}
        mv {params.fwd} {input.fwd}
        touch {output.msg}
        """

rule paladin_prepare:
    input:
        fasta_path = config['fasta_path'],
        paladin_ref_dir = config['paladin_ref_dir']
    output:
        msg = "paladin_prepare_done.txt"
    conda:
        "config/envs/paladin_env.yml"
    shell:
        """
        bash ./scripts/prepare.sh {input.fasta_path} {input.paladin_ref_dir}
        touch {output.msg}
        """

rule paladin_align:
    input:
        paladin_ref = config["paladin_db"],
        fastq_path = "{output}/{study}/{sample}/{sample}_1.fastq.gz",
        qc_complete =  "{output}/{study}/{sample}/qc_done.txt",
        prep_complete = "paladin_prepare_done.txt"
    params:
        paladin_out = "{output}/{study}/{sample}/{sample}"
    output:
        sam_out = "{output}/{study}/{sample}/{sample}.sam",
        txt = "{output}/{study}/{sample}/paladin_done.txt"
    conda:
        "config/envs/paladin_env.yml"
    resources:
        ncores = cpus_paladin
    shell:
        """
        bash ./scripts/align.sh {input.paladin_ref} {input.fastq_path} {params.paladin_out} {resources.ncores}
        touch {output.txt}
        rm {input.fastq_path}
        """
