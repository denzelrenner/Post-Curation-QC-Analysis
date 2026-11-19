import os

rule find_gap_positions:
    input:
        fa='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/corrected_fasta/{sample}_Haplome{haplomeNumber}.fa',
    output:
        gaps='{sample}/s7_overlay_gaps/{sample}_Haplome{haplomeNumber}/{sample}_Haplome{haplomeNumber}_gap_positions.tsv',
    resources:
        runtime=60,
        mem_mb=2000,
    threads:
        2
    localrule:
        True
    wildcard_constraints:
        haplomeNumber=r"\d+"
    conda:
        "../../envs/qc-tools.yml"
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # run gfastats
        gfastats {input.fa} -b gaps -j {threads} > {output.gaps}
        '''


