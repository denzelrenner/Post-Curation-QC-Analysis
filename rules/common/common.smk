# helper function to be used in the actual snakefile, and in other snakemake files as well
import os
from datetime import datetime

# wildcard constaints
wildcard_constraints:
    sample='|'.join(SAMPLES), # limit value of sample to our input bams


### GENERAL ###

# get current date which will be used to name combined QC files
CURRENTDATE=datetime.today().strftime('%d_%m_%Y')


def set_rule_target():

    outFiles = []

    # get target files for basic stats qc
    for smp in SAMPLES:
        outFiles.append(f'{smp}/s3_basic_stats_qc/qc_summary/{smp}_basicstats_qc.xlsx')

    # get target files for compleasm dbs
    outFiles.extend(expand(f"{config['compleasm_download_path']}/{{lineage}}.done",lineage=ALLCOMPLEASMLINEAGES))

    return outFiles


