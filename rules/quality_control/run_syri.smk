import itertools
def get_syri_comparisons(w):

    # set out command var
    outCommands = []
    
    # get input population fastas
    inputFastas = expand(f'{w.sample}/s2_rename_reorient_fastas/{w.sample}_Haplome{{hap}}/corrected_fasta/{w.sample}_Haplome{{hap}}.fa',hap=sampleInfo[w.sample]['ploidyAsArray'])

    # add _ to everything but the last element
    inputFastas = [f'{inputFastas[i]}' for i in range(len(inputFastas))]

    # add input fasta to comparison list
    comparisonFiles = [sampleInfo[w.sample]['reference_fasta']]  + inputFastas

    # create pairwise file comparisions
    pairwise_comparisons = itertools.pairwise(comparisonFiles)

    return [f'{ref}_{query}' for tup in pairwise_comparisons for ref,query in tup]

def create_genomes_txt(w):

    # out text var
    outText = ''
    
    # get input population fastas
    inputFastas = expand(f'{w.sample}/s2_rename_reorient_fastas/{w.sample}_Haplome{{hap}}/corrected_fasta/{w.sample}_Haplome{{hap}}.fa',hap=sampleInfo[w.sample]['ploidyAsArray'])

    # combine lists
    comparisonFiles = [sampleInfo[w.sample]['reference_fasta']]  + inputFastas

    # loop through everything and add their string that will go in the file
    for genomeFasta in comparisonFiles:

        # get name of the haplome
        genomeName = genomeFasta.split('/')[-1].removesuffix('.fa')

        # create its input line
        lineInTxt = f'{genomeFasta}\t{genomeName}\tlw:1.5\n'

        # add to final string
        outText += lineInTxt

    return outText

def get_syri_plotsr_sr_parameter(w):

    # variable with params
    srParameter = "" 

    # get comparisons
    comparisons = get_syri_comparisons(w)

    # loop through each comparison and build the path to the .out file out of syric
    for c in comparisons:

        # get individual ref, query fastas
        refFa,queryFa = c.split('_')

        # get name for haplome
        refName,queryName = refFa.split('/')[-1].removesuffix('.fa'),queryFa.split('/')[-1].removesuffix('.fa')

        # get the comparison prefix 
        comparison = f"{refName}_{queryName}"

        # get the out file
        outFile = f"--sr {w.sample}/s6_run_syri/run_syri/{comparison}/{comparison}syri.out "

        # add to parameter var
        srParameter += outFile

    return srParameter

rule pairwise_align_for_syri:
    input:
        fa=lambda w: expand(f'{w.sample}/s2_rename_reorient_fastas/{w.sample}_Haplome{{hap}}/corrected_fasta/{w.sample}_Haplome{{hap}}.fa',hap=sampleInfo[w.sample]['ploidyAsArray']),
    output:
        done='{sample}/s6_run_syri/pairwise_alignments/.done.txt',
    resources:
        runtime=480,
        mem_mb=150000
    threads:
        int(workflow.cores * 0.5)
    params:
        comparisions=get_syri_comparisons,
        mainOutdir='{sample}/s6_run_syri/pairwise_alignments'
    conda:
        "../../envs/mapping-tools.yml"
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # if outdir exists get rid of it
        if [ -d {params.mainOutdir} ];then
            rm -r {params.mainOutdir}
            mkdir -p {params.mainOutdir}
        fi

        # loop through comparisons
        for c in {params.comparisons};do

            # assign hap 1 and 2 to variables by splitting input string
            IFS='_' read -r ref_fa query_fa <<< $c

            ref_prefix=$(basename ${{hap1_raw}} .fa)
            query_prefix=$(basename ${{hap2_raw}} .fa)

            # make outdir
            comparison="${{ref_prefix}}_${{query_prefix}}"
            outdir={params.mainOutdir}/$comparison
            mkdir -p $outdir
            cd $outdir

            # run minimap command
            minimap2 -ax asm5 -t {threads} --eqx $ref_fa $query_fa | samtools sort -@ {threads} -O BAM - > "${{comparison}}.bam"
            
            # index bam
            samtools index "${{comparison}}.bam"

        done

        # move into main outdir
        cd {params.mainOutdir}

        echo 'done' > {output}
        '''

rule run_syri:
    input:
        '{sample}/s6_run_syri/pairwise_alignments/.done.txt',
    output:
        '{sample}/s6_run_syri/run_syri/.done.txt',
    resources:
        runtime=480,
        mem_mb=150000,
    params:
        comparisions=get_syri_comparisons,
        mainOutdir='{sample}/s6_run_syri/run_syri',
    threads:
        int(workflow.cores * 0.5)
    conda:
        "../../envs/syri.yml"
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # if outdir exists get rid of it
        if [ -d {params.mainOutdir} ];then
            rm -r {params.mainOutdir}
            mkdir -p {params.mainOutdir}
        fi

        # loop through comparisons
        for c in {params.comparisons};do

            # assign hap 1 and 2 to variables by splitting input string
            IFS='_' read -r ref_fa query_fa <<< $c

            ref_prefix=$(basename ${{hap1_raw}} .fa)
            query_prefix=$(basename ${{hap2_raw}} .fa)

            # make outdir
            comparison="${{ref_prefix}}_${{query_prefix}}"
            outdir={params.mainOutdir}/$comparison
            mkdir -p $outdir
            cd $outdir

            # run syri
            syri -c $outdir/"${{comparison}}.bam" --nc {threads} -r $ref_fa -q $query_fa -k -F B --prefix $comparison

        done

        # move into main outdir
        cd {params.mainOutdir}

        echo 'done' > {output}    	
        '''

rule run_plotsr:
    input:
        '{sample}/s6_run_syri/run_syri/.done.txt',
    output:
        '{sample}/s6_run_syri/generate_synteny_plot/{sample}_ref_haps_synteny_plot.pdf',
    resources:
        runtime=480,
        mem_mb=150000
    threads:
        int(workflow.cores * 0.5)
    conda:
        "../../envs/plotsr.yml"
    params:
        mainOutdir='{sample}/s6_run_syri/generate_synteny_plot',
        genomeTxt=create_genomes_txt,
        srParam=get_syri_plotsr_sr_parameter,
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # if outdir exists get rid of it
        if [ -d {params.mainOutdir} ];then
            rm -r {params.mainOutdir}
            mkdir -p {params.mainOutdir}
        fi

        # move into final outdir
        cd {params.mainOutdir}

        cat << EOF > genomes.txt
{params.genomeTxt}
EOF

        plotsr {params.srParam}--genomes genomes.txt -H 8 -W 5 -o {wildcards.sample}_ref_haps_synteny_plot.pdf
        '''