
OUTDIR=config['outdir']
sample=config['sample']
threads=config['threads']
dropseq='tools/Drop-seq_tools-2.5.3'
picard='tools/picard.jar'
bcl2fastq='tools/configureBclToFastq.pl'

rule all:
    input:
        human=unpack('{OUTDIR}/star/{sample}_human_final.bam')
        mouse=unpack('{OUTDIR}/star/{sample}_mouse_final.bam')
    threads:
        config['threads']

rule FastqtoSam:
    input:
        read1=config['r1']
        read2=config['r2']
    output:
        '{OUTDIR}/{sample}/mapping/unaligned.bam'
    params:
        "PLATFORM=illumina SORT_ORDER=queryname SAMPLE_NAME={sample}"
    shell:
        """
        java -jar {picard} FastqtoSam F1={input.read2} F2={input.read1} O={output} {params}
        """

rule tagCellBarcodes:
    input:
        '{OUTDIR}/{sample}/mapping/unaligned.bam
    output:
        '{OUTDIR}/{sample}/mapping/unaligned_tagged_cell.bam',
        summary='{OUTDIR}/{sample}/mapping/unaligned_tagged_cell.summary'
    params:
        "BASE_RANGE=0-16 BASE_QUALITY=10 BARCODE_READ=1 DISCARD_READ=False TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1"
    shell:
        """
        {dropseq}/TagBamWithReadSequenceExtended INPUT={input.human} OUTPUT={output.human} \
        SUMMARY={output.human_summary} {params}

        {dropseq}/TagBamWithReadSequenceExtended INPUT={input.mouse} OUTPUT={output.mouse} \
        SUMMARY={output.mouse_summary} {params}
        """
rule tagUMI:
    input:
        human='{OUTDIR}/{sample}/mapping/unaligned_tagged_cell_human.bam'
        mouse='{OUTDIR}/{sample}/mapping/unaligned_tagged_cell_mouse.bam'
    output:
        human=temp('{OUTDIR}/{sample}/mapping/unaligned_tagged_molecular_human.bam')
        human_summary=
        mouse=temp('{OUTDIR}/{sample}/mapping/unaligned_tagged_molecular_mouse.bam')
        mouse_summary=
    params:
        "BASE_RANGE=16-28 BASE_QUALITY=10 BARCODE_READ=1 DISCARD_READ=False TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1"
    shell:
        """
        {dropseq}/TagBamWithReadSequenceExtended INPUT={input.human} OUTPUT={output.human} \
        SUMMARY={output.human_summary} {params}

        {dropseq}/TagBamWithReadSequenceExtended INPUT={input.mouse} OUTPUT={output.mouse} \
        SUMMARY={output.mouse_summary} {params}
        """

rule TrimStartingSequence:
    input:
        human='{OUTDIR}/{sample}/mapping/unaligned_tagged_molecular_human.bam'
        mouse='{OUTDIR}/{sample}/mapping/unaligned_tagged_molecular_mouse.bam'
    output:
        human=temp('{OUTDIR}/{sample}/mapping/tagged_trimmed_human.bam')
        mouse=temp('{OUTDIR}/{sample}/mapping/tagged_trimmed_mouse.bam')
    params:
        "SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5"
    shell:
        """
        {dropseq}/TrimStartingSequence INPUT={input.human} OUTPUT={output.human} OUTPUT_SUMMARY={human.summary} {params}
        {dropseq}/TrimStartingSequence INPUT={input.human} OUTPUT={output.human} OUTPUT_SUMMARY={human.summary} {params}
        """

rule PolyATrimmer:
    input:
        human=temp('{OUTDIR}/{sample}/mapping/tagged_trimmed_human.bam')
        mouse=temp('{OUTDIR}/{sample}/mapping/tagged_trimmed_mouse.bam')
    output:
        human=temp('{OUTDIR}/{sample}/mapping/human_tagged_polyA_adapter_trimmed.bam')
        mouse=temp('{OUTDIR}/{sample}/mapping/mouse_tagged_polyA_adapter_trimmed.bam')
    params:
        "MISMATCHES=0 NUM_BASES=6"
    shell:
        """
        {dropseq}/PolyATrimmer OUTPUT_SUMMARY={reports_dir}/human_remove_polyA_report.txt INPUT={input.human} OUTPUT={output.human} {params}
        {dropseq}/PolyATrimmer OUTPUT_SUMMARY={reports_dir}/mouse_remove_polyA_report.txt INPUT={input.mouse} OUTPUT={output.mouse} {params}
        """


rule SamtoFastq:
    input:
        human='{OUTDIR}/{sample}/mapping/human_tagged_polyA_adapter_trimmed.bam'
        mouse='{OUTDIR}/{sample}/mapping/mouse_tagged_polyA_adapter_trimmed.bam'
    output:


rule index_genome:
    input:
        mouse_ref = config['mouse_ref']
        mouse_annotation=config['mouse_annotation']
        human_annotation=config['human_annotation']
        human_ref = config['human_ref']
    output:
        human_index=directory('species_data/human/star_index')
        mouse_index=directory('species_data/mouse/star_index')
    params:
        human_index=config['human_star_index']
        mouse_index=config['mouse_star_index']
    run:
        if snakemake.params.human_index == None:
            shell("STAR --runMode genomeGenerate \
                 --runThreadN {threads} \
                 --genomeDir {output.human_index} \
                 --genomeFastaFiles {input.human_ref} \
                 --sjdbGTFfile {input.human_annotation}")

        if snakemake.params.mouse_index == None:
            shell("STAR --runMode genomeGenerate \
                 --runThreadN {threads} \
                 --genomeDir {output.mouse_index} \
                 --genomeFastaFiles {input.mouse_ref} \
                 --sjdbGTFfile {input.mouse_annotation}")
        else:
            shell("ln -s {params.human_index} {output.human_index}")
            shell("ln -s {params.mouse_index} {output.mouse_index}")

rule STAR_align:
    input:
        human=unpack('{OUTDIR}/{sample}/mapping/human_{read}.fastq',read=['r1','r2'])
        mouse=unpack('{OUTDIR}/{sample}/mapping/mouse_{read}.fastq',read=['r1','r2'])

    output:
        human=temp('{OUTDIR}/{sample}/mapping/human_Aligned_out.bam')
        mouse=temp('{OUTDIR}/{sample}/mapping/mouse_Aligned_out.bam')
    params:
        human='species_data/human/star_index'
        mouse='species/mouse/star_index'
        human_annotation=config['human_annotation']
        mouse_annotation=config['mouse_annotation']
        aln="--readFilesCommand zcat --outSAMtype BAM Unsorted"
    threads:
        config['threads']
    shell:
        """
        STAR --genomeDir {params.human} \
            --readFilesIn {input.human} \
            --outFileNamePrefix {OUTDIR}/{sample}/mapping/human_ \
            --sjdbGTFfile {params.human_annotation} \
            {params.aln} \
            --runThreadN {threads}
        STAR --genomeDir {params.mouse} \
            --readFilesIn {input.mouse} \
            --outFileNamePrefix {OUTDIR}/{sample}/mapping/mouse_ \
            --sjdbGTFfile {params.mouse_annotation} \
            {params.aln} \
            --runThreadsN {threads}
        """


rule SortSam:
    input:
        human='{OUTDIR}/{sample}/mapping/human_Aligned_out.bam'
        mouse='{OUTDIR}/{sample}/mapping/mouse_Aligned_out.bam'
    output:
        human=temp('{OUTDIR}/{sample}/mapping/human_sorted_aligned.bam')
        mouse=temp('{OUTDIR}/{sample}/mapping/mouse_sorted_aligned.bam')
    shell:
        """
        java -jar {picard} SortSam I={input.human} O={output.human} SORT_ORDER=queryname
        java -jar {picard} SortSam I={input.mouse} O={output.mouse} SORT_ORDER=queryname
        """

rule MergeAlignment:
    input:
        human_aln='{OUTDIR}/{sample}/mapping/human_sorted_aligned.bam'
        human_ubam='{OUTDIR}/{sample}/mapping/human_tagged_polyA_adapter_trimmed.bam'
        mouse_aln='{OUTDIR}/{sample}/mapping/mouse_sorted_aligned.bam'
        mouse_ubam='{OUTDIR}/{sample}/mapping/mouse_tagged_polyA_adapter_trimmed.bam'
    output:
        human=temp('{OUTDIR}/{sample}/mapping/human_merged.bam')
        mouse=temp('{OUTDIR}/{sample}/mapping/mouse_merged.bam')
    params:
        human_ref=config['human_ref']
        mouse_ref=config['mouse_ref']
    shell:
        """
        java -jar {picard} MergeBamAlignment REFERENCE_SEQUENCE={params.human_ref} UNMAPPED_BAM={input.human_ubam} ALIGNED_BAM={input.human_aln} OUTPUT={output.human} {params}
        java -jar {picard} MergeBamAlignment REFERENCE_SEQUENCE={params.mouse_ref} UNMAPPED_BAM={input.mouse_ubam} ALIGNED_BAM={input.mouse_aln} OUTPUT={output.mouse} {params}
        """

rule TagReadWithGeneExon:
    input:
        human='{OUTDIR}/{sample}/mapping/human_merged.bam'
        mouse='{OUTDIR}/{sample}/mapping/mouse_merged.bam'
    output:
        human=temp('{OUTDIR}/{sample}/mapping/human_merged_tagged.bam')
        mouse=temp('{OUTDIR}/{sample}/mapping/mouse_merged_tagged.bam')
    shell:
        """
        {dropseq}/TagReadWithGeneExon I={input.human} O={output.human} ANNOTATIONS_FILE={params.human_annotaiton} TAG=GE
        {dropseq}/TagReadWithGeneExon I={input.human} O={output.human} ANNOTATIONS_FILE={params.human_annotaiton} TAG=GE
        """

rule filer_mm_reads:
    input:
        human='{OUTDIR}/{sample}/mapping/human_merged_tagged.bam'
        mouse='{OUTDIR}/{sample}/mapping/mouse_merged_tagged.bam'
    output:
        human='{OUTDIR}/{sample}/mapping/human_final.bam'
        mouse='{OUTDIR}/{sample}/mapping/mouse_final.bam'
    shell:
        """
        python scripts/filter_mm_reads.py --in-bam {input.human} --out-bam {output.human}
        python scripts/filter_mm_reads.py --in-bam {input.mouse} --out-bam {output.mouse}
        """
