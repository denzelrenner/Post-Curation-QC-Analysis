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

    # here we cant use expand because each sample doesnt have the same number of haplotypes
    for smp in SAMPLES:

        # mapped reads target file
        outFiles.extend(expand(f"{smp}/s5_map_reads_to_assem/{smp}_Haplome{{hap}}/hifi/{smp}_Haplome{{hap}}_{{suff}}",hap=sampleInfo[smp]['ploidyAsArray'],suff=['hifi_sorted.bam','hifi_sorted.bam.bai','hifi_to_hifiasm.tsv']))

        # gap positions target file
        outFiles.extend(expand(f'{smp}/s7_overlay_gaps/{smp}_Haplome{{hap}}/{smp}_Haplome{{hap}}_gap_positions.tsv',hap=sampleInfo[smp]['ploidyAsArray']))


    # get target files for compleasm dbs
    outFiles.extend(expand(f"{config['compleasm_download_path']}/{{lineage}}.done",lineage=ALLCOMPLEASMLINEAGES))

    # get target files for basic stats qc
    outFiles.extend(expand('{smp}/s3_basic_stats_qc/qc_summary/{smp}_basicstats_qc.xlsx',smp=SAMPLES))

    # get target file for syri synteny plot
    outFiles.extend(expand('{smp}/s6_run_syri/generate_synteny_plot/{smp}_ref_haps_synteny_plot.pdf',smp=SAMPLES))

    return outFiles


