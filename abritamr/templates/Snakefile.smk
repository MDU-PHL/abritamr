import pathlib
configfile: 'config.yaml'
workdir: config['workdir']


SAMPLES = config['samples'].split()
SINGULARITY_PATH = config['singularity_path']
FINAL_OUTPUT = config['final_output']
MDU_QC = config['mduqc']
SCRIPT_PATH = config['script_path']
QC = config['qc']
SPECIES_DETECT = config['species_detect'] == 'species'
RUN_AMRFINDER = config['amrfinder'] == 'run_amrfinder'
DB = config['database_version']

def get_contigs(wildcards):
    return {"contigs": f"{pathlib.Path(wildcards.sample, 'contigs.fa')}"}


def get_amr(wildcards):
    return {"amr": pathlib.Path(wildcards.sample, f"{wildcards.sample}.out").absolute().as_posix()}


rule all:
    input: FINAL_OUTPUT

if RUN_AMRFINDER:
    rule amr_finder:
        input: unpack(get_contigs)
        output: "{sample}/{sample}.out"
        singularity: SINGULARITY_PATH
        conda: "amrfinderplus.yaml"
        shell: "amrfinder -n {input} -o {output} -t 1"

rule collate:
    input: expand("{sample}/{sample}.out", sample=SAMPLES)
    output: FINAL_OUTPUT
    params:
        mduqc = MDU_QC,
        script_path = SCRIPT_PATH,
        qc = QC,
        db = DB
    shell:
        """
        python3 "{params.script_path}/collate.py" {params.db} {params.mduqc} {params.qc} 
        """
