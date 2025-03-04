import os

configfile: 'config/config.yml'
configfile: 'config/cluster.json'

cpus_ena = config['ena_download']['cores']
cpus_qc = config['metagen_qc']['cores']
cpus_paladin = config['paladin']['cores']
host_ref = config['host_ref']


INPUT_FILE = config['input_file']
OUTPUT_FASTQ_DIR = config['output_fastq_dir']
FASTA_DIR = config['fasta_dir']
FASTA = config['fasta']
OUTPUT_ALIGN_DIR = config['output_align_dir']
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
        fwd = "{output}/{study}/{sample}/{sample}_1.fastq.gz",
        rev = "{output}/{study}/{sample}/{sample}_2.fastq.gz"
    params:
        outdir = "{output}/{study}/{sample}"
    conda:
        "config/envs/ena_download.yml"
    resources:
        ncores = cpus_ena
    shell:
        """
        fastq-dl --cpus {resources.ncores} -a {wildcards.sample} --outdir {params.outdir}
        """

rule metagen_qc:
    input:
        fwd = "{output}/{study}/{sample}/{sample}_1.fastq.gz",
        rev = "{output}/{study}/{sample}/{sample}_2.fastq.gz"
    output:
        "{output}/{study}/{sample}/qc_done.txt"
    params:
        bwa_ref = host_ref,
        fwd = "{output}/{study}/{sample}/{sample}_clean_1.fastq.gz",
        rev = "{output}/{study}/{sample}/{sample}_clean_2.fastq.gz"
    conda:
        "config/envs/metagen_qc.yml"
    resources:
        ncores = cpus_qc
    shell:
        """
        ./scripts/metagen-fastqc.sh -t {resources.ncores} -c {params.bwa_ref} -f {input.fwd} -r {input.rev}
        mv {params.fwd} {input.fwd}; mv {params.rev} {input.rev}
        touch {output}
        """

rule paladin_align:
    input:
        fastq_path = "{output}/{study}/{sample}/{sample}_1.fastq.gz",
        rev = "{output}/{study}/{sample}/{sample}_2.fastq.gz",
        qc_flag =  "{output}/{study}/{sample}/qc_done.txt"
    output:
        tsv = "{output}/{study}/{sample}/paladin_done.txt"
    conda:
        "config/envs/paladin_env.yml"
    resources:
        ncores = cpus_paladin
    shell:
        """
        bash ./scripts/align.sh {FASTA_DIR} {FASTA} {input.fastq_path} {resources.ncores} {OUTPUT_ALIGN_DIR} {output.tsv}
        rm {input.fastq_path}
        rm {input.rev}
        """