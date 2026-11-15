rule map_hifi_to_assem:
    input:
        fa='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/corrected_fasta/{sample}_Haplome{haplomeNumber}.fa',
        fai='{sample}/s2_rename_reorient_fastas/{sample}_Haplome{haplomeNumber}/corrected_fasta/{sample}_Haplome{haplomeNumber}.fa.fai',
        hifi="{sample}/s4_get_input_reads/Hifi/{sample}.fastq.gz"
    output:
        bam='{sample}/s5_map_reads_to_assem/{sample}_Haplome{haplomeNumber}/hifi/{sample}_Haplome{haplomeNumber}_hifi_sorted.bam',
        bam_idx='{sample}/s5_map_reads_to_assem/{sample}_Haplome{haplomeNumber}/hifi/{sample}_Haplome{haplomeNumber}_hifi_sorted.bam.bai',
        summary_tsv='{sample}/s5_map_reads_to_assem/{sample}_Haplome{haplomeNumber}/hifi/{sample}_Haplome{haplomeNumber}_hifi_to_hifiasm.tsv',
    resources:
        runtime=480,
        mem_mb=150000
    threads:
        int(workflow.cores * 0.5)
    wildcard_constraints:
        haplomeNumber=r"\d+"
    conda:
        "../../envs/mapping-tools.yml"
    shell:
        '''
        export PATH=$CONDA_PREFIX/bin:$PATH

        # map hifi reads to fasta
        minimap2 -ax map-pb -t {threads} {input.fa} {input.hifi} | samtools sort -@{threads} -O BAM -o {output.bam} -

		# index bam file
        samtools index {output.bam} --threads {threads} -o {output.bam_idx}
        
        # run samtools flagstat to get mapping metrics
        samtools flagstat --output-fmt=tsv {output.bam} > {output.summary_tsv}
        '''



