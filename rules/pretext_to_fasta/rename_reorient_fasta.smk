rule create_uniform_name_for_input_fasta:
    input:
        done='{sample}/s1_agp_to_fasta/.done.txt',
    output:
        fa='{sample}/s1_agp_to_fasta/renamed_main_scaffold_file/{sample}.hap{haplomeNumber}.fa',
    params:
        ploidy=lambda w: sampleInfo[w.sample]['ploidy']
    threads:
        1
    resources:
        mem_mb=1000,
        runtime=20,
    wildcard_constraints:
        haplomeNumber=r"\d+"
    localrule:
        True
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # depending on the ploidy the output fastas look a little different
        if [ "{params.ploidy}" -eq "2" ];then

            mv {wildcards.sample}/s1_agp_to_fasta/{wildcards.sample}.hap{wildcards.haplomeNumber}.1*.fa {output.fa}

        elif [ "{params.ploidy}" -gt "2" ];then

            mv {wildcards.sample}/s1_agp_to_fasta/{wildcards.sample}.1.hap{wildcards.haplomeNumber}s*.fa {output.fa}

        fi
        '''

rule align_fasta_to_ref:
    input:
        fa='{sample}/s1_agp_to_fasta/renamed_main_scaffold_file/{sample}.hap{haplomeNumber}.fa',
        done='{sample}/s1_agp_to_fasta/.done.txt',
    output:
        delta='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/align_fasta_to_ref/output.delta',
    conda:
        "../../envs/mummer.yml"
    params:
        ref=lambda w: sampleInfo[w.sample]['reference_fasta'],
    threads:
        8
    resources:
        mem_mb=20000,
        runtime=720,
    wildcard_constraints:
        haplomeNumber=r"\d+"
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # align fasta to ref
        nucmer -t {threads} --prefix={wildcards.sample}/s2_rename_reorient_fastas/{wildcards.sample}_Haplome{wildcards.haplomeNumber}/align_fasta_to_ref/output {params.ref} {input.fa}
        '''

rule create_coords_file:
    input:
        delta='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/align_fasta_to_ref/output.delta'
    output:
        coords='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/align_fasta_to_ref/output.coords',
        coords_idx='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/align_fasta_to_ref/output.coords.idx',
        deltagz='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/align_fasta_to_ref/output.uniqueAnchorFiltered_l10000.delta.gz',
    conda:
        "../../envs/pd-tools.yml"
    params:
        dotPrepScript=os.path.join(config['script_dir'],'python_scripts','DotPrep.py')
    threads:
        2
    localrule:
        True
    resources:
        mem_mb=1000,
        runtime=120,
    wildcard_constraints:
        haplomeNumber=r"\d+"
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH
        
        # needs python 3.6
        {params.dotPrepScript} --delta {output.delta} --out {wildcards.sample}/s2_rename_reorient_fastas/{wildcards.sample}_Haplome{wildcards.haplomeNumber}/align_fasta_to_ref/output
        '''

rule rename_reorient_fasta:
    input:
        fa='{sample}/s1_agp_to_fasta/renamed_main_scaffold_file/{sample}.hap{haplomeNumber}.fa',
        coords_idx='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/align_fasta_to_ref/output.coords.idx',
    output:
        outFa='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/corrected_fasta/{sample}_Haplome{haplomeNumber}.fa',
        outXlsx='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/corrected_fasta/{sample}_Haplome{haplomeNumber}_remapping_sense_check.xlsx',
    conda:
        "../../envs/pd-tools.yml",
    params:
        ploidy=lambda w: sampleInfo[w.sample]['ploidy'],
        qcScript=os.path.join(config['script_dir'],'python_scripts','reorder_reorient_chr_to_ref.py'),
        ref=lambda w: sampleInfo[w.sample]['reference_fasta'],
        xlsxPrefix=lambda w:f'{w.sample}_Haplome{w.haplomeNumber}_remapping_sense_check',
        outdir=lambda w:f'{w.sample}/s2_rename_reorient_fastas/{w.sample}_Haplome{w.haplomeNumber}/corrected_fasta'
    resources:
        mem_mb=1000,
        runtime=30,
    localrule:
        True
    wildcard_constraints:
        haplomeNumber=r"\d+"
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        python3 {params.qcScript} -s {wildcards.sample} --sample_fasta {input.fa} --ref_fasta {params.ref} -nhap {wildcards.haplomeNumber} -c {input.coords_idx} -o {params.xlsxPrefix} -od {params.outdir}
        '''
