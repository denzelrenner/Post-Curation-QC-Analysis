import os

rule combine_debris_fasta:
    input:
        fa=lambda w:expand(f'{w.sample}/s1_agp_to_fasta/renamed_main_scaffold_file/{w.sample}.hap{{hapNum}}.fa',hapNum=sampleInfo[w.sample]['ploidyAsArray']),
        step1Done='{sample}/s1_agp_to_fasta/.done.txt', 
    output:
        debrisFasta='{sample}/s2_rename_reorient_fastas/{sample}_HaplomeDebris/corrected_fasta/{sample}_HaplomeDebris.fa',
        senseCheck='{sample}/s2_rename_reorient_fastas/{sample}_HaplomeDebris/corrected_fasta/.debris_fasta_files.txt',
    resources:
        mem_mb=1000,
        runtime=30,
    threads:
        1
    params:
        debrisHapNumber=lambda w:sampleInfo[w.sample]['ploidy'] + 1,
    localrule:
        True
    shell:
        '''
        # combine debris fastas
	ls {wildcards.sample}/s1_agp_to_fasta/*.fa > {output.senseCheck}
        cat {wildcards.sample}/s1_agp_to_fasta/*.fa > {output.debrisFasta}
        '''
