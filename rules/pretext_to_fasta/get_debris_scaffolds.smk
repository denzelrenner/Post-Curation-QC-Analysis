import os
def get_scaffolding_debris(w):

    if sampleInfo[w.sample]['ploidy'] == 2:

        # get main scaffolds
        mainFastas = expand(f'{w.sample}/s1_agp_to_fasta/{w.sample}.hap{{hap}}.1.primary.curated.fa',hap=sampleInfo[w.sample]['ploidyAsArray'])

        # get all fastas
        allFastas = [f for f in os.listdir(f'{w.sample}/s1_agp_to_fasta') if f.endswith('.curated.fa')]

        # get debris only
        debrisFastas = set(allFastas).difference(set(mainFastas))

        return list(debrisFastas)

    elif sampleInfo[w.sample]['ploidy'] > 2:

        # get main scaffolds
        mainFastas = expand(f'{w.sample}/s1_agp_to_fasta/{w.sample}1.hap{{hap}}s.curated.fa',hap=sampleInfo[w.sample]['ploidyAsArray'])

        # get all fastas
        allFastas = [f for f in os.listdir(f'{w.sample}/s1_agp_to_fasta') if f.endswith('.curated.fa')]

        # get debris only
        debrisFastas = set(allFastas).difference(set(mainFastas))

        return list(debrisFastas)


rule combine_debris_fasta:
    input:
        fa=get_scaffolding_debris,
        step1Done='{sample}/s1_agp_to_fasta/.done.txt', 
    output:
        '{sample}/s2_rename_reorient_fastas/{sample}_HaplomeDebris/corrected_fasta/{sample}_HaplomeDebris.fa',
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
        # make debris haplome dir
        mkdir -p {wildcards.sample}/s2_rename_reorient_fastas/{wildcards.sample}_Haplome{params.debrisHapNumber}/corrected_fasta

        # combine debris fastas
        cat {input.fa} > {output}
        '''
