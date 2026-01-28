import os

configfile: 'config/config.yml'
configfile: 'config/cluster.json'

cpus_ena = config['ena_download']['cores']
cpus_qc = config['metagen_qc']['cores']
cpus_paladin = config['paladin']['cores']
host_ref = config['host_ref']

INPUT_FILE = config['input_file']
OUTPUT_FASTQ_DIR = config['output_fastq_dir']
OUPUT_SUMMARY_DIR = config['output_summary_dir']
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


logs_dirname = "logs/" + config['logfile_name']
if not os.path.exists(logs_dirname):
    os.makedirs(logs_dirname) 

rule targets:
   input:
       expand([OUPUT_SUMMARY_DIR+"{study}/{sample}_tsv_to_summary_done.txt"], zip, study=STUDIES, sample=SAMPLES)

rule ena_download:
    output:
        msg = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "ena_done.txt")
    params:
        fwd = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "{sample}_1.fastq.gz"), 
        outdir = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}"),
        sample = "{sample}",
        log_file = logs_dirname + "/ena_download.{sample}.err"
    conda:
        "config/envs/ena_download.yml"
    resources:
        ncores = cpus_ena, 
        load = 2000
    shell:
        """
        bash scripts/ena_download_wrapper.sh \
            {params.sample} \
            {params.outdir} \
            {resources.ncores} \
            {params.fwd}
        """

rule metagen_qc:
    input:
        ena_done = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "ena_done.txt")
    output:
        msg = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "qc_done.txt")
    params:
        bwa_ref = host_ref,
        fwd = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "{sample}_1.fastq.gz"), 
        fwd_clean = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "{sample}_1_clean.fastq.gz"),
        log_file = os.path.join(logs_dirname, "metagen_qc.{sample}.err")
    conda:
        "config/envs/metagen_qc.yml"
    resources:
        ncores = cpus_qc
    shell:
        """
        ./scripts/metagen-fastqc.sh -t {resources.ncores} -c {params.bwa_ref} -f {params.fwd}
        mv {params.fwd_clean} {params.fwd}
        touch {output.msg}
        """

rule paladin_prepare:
    input:
        fasta_path = config['fasta_path'],
        paladin_ref_dir = config['paladin_ref_dir']
    output:
        msg = os.path.join(PALADIN_REF_DIR, "paladin_prepare_done.txt"),
        index = os.path.join(PALADIN_REF_DIR, FASTA_NAME)
    log:
        log_file = os.path.join(logs_dirname, "paladin_prepare.err")
    conda:
        "config/envs/paladin_env.yml"
    shell:
        """
        bash ./scripts/prepare.sh {input.fasta_path} {input.paladin_ref_dir}
        touch {output.msg}
        """

rule paladin_align:
    input:
        prep_fasta_path = os.path.join(PALADIN_REF_DIR, FASTA_NAME),
        qc_complete = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "qc_done.txt"),
        prep_complete = os.path.join(PALADIN_REF_DIR, "paladin_prepare_done.txt")
    params:
        fastq_path = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "{sample}_1.fastq.gz"),
        paladin_out = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "{sample}"),
        log_file = os.path.join(logs_dirname, "paladin_align.{sample}.err")
    output:
        msg = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "paladin_done.txt")
    conda:
        "config/envs/paladin_env.yml"
    resources:
        ncores = cpus_paladin, 
        load = 4
    shell:
        """
        #bash ./scripts/align.sh {input.prep_fasta_path} {params.fastq_path} {params.paladin_out} {resources.ncores}
        bash scripts/paladin_align_wrapper.sh \
            {input.prep_fasta_path} \
            {params.fastq_path} \
            {params.paladin_out} \
            {resources.ncores} \
            {output.msg}
        """

rule bam_to_tsv:
    input:
        paladin_complete = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "paladin_done.txt")
    output:
        msg = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "bam_to_tsv_done.txt")
    params:
        bam = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "{sample}.bam"), 
        output_dir = os.path.join(OUTPUT_FASTQ_DIR, "tsv_outputs", "{study}"), 
        log_file = os.path.join(logs_dirname, "bam_to_tsv.{sample}.err")
    conda:
        "config/envs/r_env.yml"
    shell:
        """
        mkdir -p {params.output_dir}
        bash ./scripts/bam_to_tsv_wrapper.sh {params.bam} {params.output_dir}
        touch {output.msg}
        """

rule tsv_to_summary: 
    input:
        msg = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "bam_to_tsv_done.txt")
    output:
        tsv = os.path.join(OUPUT_SUMMARY_DIR, "{study}", "{sample}_summary.tsv"), 
        log_out = os.path.join(OUPUT_SUMMARY_DIR, "{study}", "{sample}_log.tsv"), 
        msg = os.path.join(OUPUT_SUMMARY_DIR, "{study}", "{sample}_tsv_to_summary_done.txt")
    params:
        bam = os.path.join(OUTPUT_FASTQ_DIR, "{study}", "{sample}", "{sample}.bam"), 
        tsv = os.path.join(OUTPUT_FASTQ_DIR, "tsv_outputs", "{study}", "{sample}.tsv"), 
        output_dir = os.path.join(OUPUT_SUMMARY_DIR, "{study}"), 
        log_file = os.path.join(logs_dirname, "tsv_to_summary.err")
    conda:
        "config/envs/r_env.yml"
    shell:
        """
        mkdir -p {params.output_dir}
        
        # Run the wrapper script
        bash ./scripts/tsv_to_summary_wrapper.sh {params.tsv} {params.output_dir}
        
        if [[ -f "{output.tsv}" && -f "{output.log_out}" ]]; then
            rm {params.tsv}
            rm {params.bam}
            touch {output.msg}
        fi
        """