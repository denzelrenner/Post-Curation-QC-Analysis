rule index_scaffold_fasta:
    input:
        fa=lambda w: sampleInfo[w.sample]['scaffold_fasta'],
    output:
        fai='{sample}/s1_agp_to_fasta/.indexed.txt'
    conda:
        "../../envs/mapping-tools.yml"
    resources:
        mem_mb=1000,
        runtime=30,
    threads:
        4
    localrule:
        True
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # index scaffold fasta
        samtools faidx {input.fa}

        echo 'done' > {output.fai}
        '''


rule agp_to_fasta:
    input:
        agp=lambda w: sampleInfo[w.sample]['agp'],
        fa=lambda w: sampleInfo[w.sample]['scaffold_fasta'],
        indexDone='{sample}/s1_agp_to_fasta/.indexed.txt'
    output:
        done='{sample}/s1_agp_to_fasta/.done.txt'
    conda:
        "../../envs/agp-utils.yml"
    params:
        ploidy=lambda w: sampleInfo[w.sample]['ploidy'],
        pretextAsmScript=os.path.join(config['script_dir'],'python_scripts','DotPrep.py'),
        outdir='{sample}/s1_agp_to_fasta',
    log:
        out='{sample}/s1_agp_to_fasta/pretext_to_asm.log',
    resources:
        mem_mb=1000,
        runtime=30,
    localrule:
        True
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # pretext to asm
        pretext-to-asm -a {input.fa} -p {input.agp} -c Chr -o {params.outdir}/{wildcards.sample}.1.fa > {log.out} 2>&1

        # depending on the ploidy the output fastas look a little different
        if [ "{params.ploidy}" -eq "2" ];then

            fastaFiles=({params.outdir}/{wildcards.sample}.hap*.1.*.fa)
            numOfFastaFiles=${{#fastaFiles[@]}}

            if [ "$numOfFastaFiles" -eq "{params.ploidy}" ];then
                echo 'All Fasta Files Present' > {output}
            fi

        elif [ "{params.ploidy}" -gt "2" ];then

            fastaFiles=({params.outdir}/{wildcards.sample}.1.hap*s.*.fa)
            numOfFastaFiles=${{#fastaFiles[@]}}

            if [ "$numOfFastaFiles" -eq "{params.ploidy}" ];then
                echo 'All Fasta Files Present' > {output}
            fi
        
        fi
        '''


