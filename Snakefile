import configparser
import sys
import os
import pandas as pd

# set config ini path
# path is relative and set relative to working directory
# you can set environemntal variable called CONFIG_INI_PATH if the file 'config.ini' is not in your working directory
configIniPath = os.environ.get('CONFIG_INI_PATH','config.ini')
SNAKEMAKE_OUTPUTDIR = os.environ.get('SNAKEMAKE_OUTPUTDIR','')
ANALYSIS_NAME = os.environ.get('ANALYSIS_NAME','analysis')

# load in config file, 
configIni = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))

# read in config gile
configIni.read(configIniPath)

# build snakemake config
### DEFAULTS ###
config['samples_file'] = configIni['DEFAULT'].get('samples_file',fallback='samples.csv')
config['rsync_folder'] = configIni['DEFAULT'].get('rsync_folder',fallback='compiled_results')

### QUALITY CONTROL
config['compleasm_download_path'] = configIni['QUALITY_CONTROL'].get('compleasm_download_path',fallback='')

### SLURM ARGS
config['max_mem_mb_per_node'] = configIni['SLURM_ARGS'].getint('max_mem_mb_per_node',fallback=32000)
config['max_cpus_per_node'] = configIni['SLURM_ARGS'].getint('max_cpus_per_node',fallback=32)
config['slurm_partition'] = configIni['SLURM_ARGS'].get('slurm_partition',fallback='pelle')
config['slurm_account'] = configIni['SLURM_ARGS'].get('slurm_account',fallback='slu2026-1-5')
config['tmpdir'] = configIni['SLURM_ARGS'].get('tmpdir',fallback='/scratch')

### set other global variables
config['script_dir'] = os.path.join(workflow.basedir,'scripts')
config['shared_data_dir'] = os.path.join(workflow.basedir,'shared_data')



# include rules

# shared rules
include:"rules/common/input_csv_processing.smk"
include:"rules/common/common.smk"
#include:"rules/common/collect_results.smk"

### PREPARING INPUT DATA ###
# bam to fastq
include:"rules/prepare_input_data/get_hifi.smk"

# merge HiC
include:"rules/prepare_input_data/get_hic.smk"

### PRETEXT TO FASTA ###
include:"rules/pretext_to_fasta/get_debris_scaffolds.smk"
include:"rules/pretext_to_fasta/pretext_to_asm.smk"
include:"rules/pretext_to_fasta/rename_reorient_fasta.smk"

### BASICSTATS QC ###
include:"rules/quality_control/run_basicstats_qc.smk"


# all rule
rule all:
    input:
        set_rule_target()
