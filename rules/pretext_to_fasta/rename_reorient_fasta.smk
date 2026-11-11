rule align_fasta_to_ref:
    input:
        fa=lambda w: f'{w.sample}/s1_agp_to_fasta/{w.sample}.hap{{haplomeNumber}}.1.primary.curated.fa' if sampleInfo[w.sample]['ploidy'] == 2 else f'{w.sample}/s1_agp_to_fasta/{w.sample}.1.hap{{haplomeNumber}}s.curated.fa',
        done='{sample}/s1_agp_to_fasta/.done.txt',
    output:
        coords='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/align_fasta_to_ref/output.coords',
        coords_idx='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/align_fasta_to_ref/output.coords.idx',
        delta='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/align_fasta_to_ref/output.delta',
        deltagz='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/align_fasta_to_ref/output.uniqueAnchorFiltered_l10000.delta.gz',
    conda:
        "../../envs/mummer.yml"
    params:
        ref=lambda w: sampleInfo[w.sample]['reference_fasta'],
    threads:
        8
    resources:
        mem_mb=20000,
        runtime=720,
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # align fasta to ref
        nucmer -t {threads} --prefix=output {params.ref} {input.fa}

        # needs python 3.6
        DotPrep.py --delta {output.delta}
        '''

rule rename_reorient_fasta:
    input:
        fa=lambda w: f'{w.sample}/s1_agp_to_fasta/{w.sample}.hap{{haplomeNumber}}.1.primary.curated.fa' if sampleInfo[w.sample]['ploidy'] == 2 else f'{w.sample}/s1_agp_to_fasta/{w.sample}.1.hap{{haplomeNumber}}s.curated.fa',
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
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        python3 {params.qcScript} -s {wildcards.sample} --sample_fasta {input.fa} --ref_fasta {params.ref} -nhap {wildcards.haplomeNumber} -c {input.coords_idx} -o {params.xlsxPrefix} -od {params.outdir}
        '''
