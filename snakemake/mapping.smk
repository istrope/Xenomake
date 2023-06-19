OUTDIR=config['outdir']
sample=config['sample']
threads=config['threads']
dropseq=config['dropseq_tools']
picard=config['picard']


rule FastqtoSam:
    input:
        read1=config['r1'],
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
        '{OUTDIR}/{sample}/mapping/unaligned.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/mapping/unaligned_tagged_cell.bam'),
        summary='{OUTDIR}/{sample}/mapping/tag_cell_barcodes.summary'
    params:
        "BASE_RANGE=0-16 BASE_QUALITY=10 BARCODE_READ=1 DISCARD_READ=False TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1"
    shell:
        """
        {dropseq}/TagBamWithReadSequenceExtended INPUT={input} OUTPUT={output.bam} \
        SUMMARY={output.summary} {params}
        """
rule tagUMI:
    input:
        '{OUTDIR}/{sample}/mapping/unaligned_tagged_cell.bam',
    output:
        bam=temp('{OUTDIR}/{sample}/mapping/unaligned_tagged_molecular.bam'),
        summary='{OUTDIR}/{sample}/mapping/tag_umi.summary',
    params:
        "BASE_RANGE=16-28 BASE_QUALITY=10 BARCODE_READ=1 DISCARD_READ=False TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1"
    shell:
        """
        {dropseq}/TagBamWithReadSequenceExtended INPUT={input} OUTPUT={output.bam} SUMMARY={output.summary} {params}
        """

rule TrimStartingSequence:
    input:
        '{OUTDIR}/{sample}/mapping/unaligned_tagged_molecular.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/mapping/tagged_trimmed.bam'),
        summary='{OUTDIR}/{sample}/mapping/trimmed_starting_sequence.summary'
    params:
        "SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5"
    shell:
        """
        {dropseq}/TrimStartingSequence INPUT={input} OUTPUT={output.bam} OUTPUT_SUMMARY={output.summary} {params}
        """

rule PolyATrimmer:
    input:
        '{OUTDIR}/{sample}/mapping/tagged_trimmed.bam'
    output:
        bam=temp('{OUTDIR}/{sample}/mapping/tagged_polyA_adapter_trimmed.bam'),
        summary='{OUTDIR}/{sample}/mapping/polyATrimmer.summary'
    params:
        "MISMATCHES=0 NUM_BASES=6"
    shell:
        """
        {dropseq}/PolyATrimmer OUTPUT_SUMMARY={output.summary} INPUT={input} OUTPUT={output.bam} {params}
        """


rule SamtoFastq:
    input:
        bam='{OUTDIR}/{sample}/mapping/tagged_polyA_adapter_trimmed.bam'
    output:
        read1='{OUTDIR}/{sample}/mapping/{sample}_processed_R1.fastq',
        read2='{OUTDIR}/{sample}/mapping/{sample}_processed_R2.fastq'
    params:
    shell:
        """
        java -jar {picard} SamToFastq I={input} FASTQ={output.read1} SECOND_END_FASTQ={output.read2}
        """


rule index_genome:
    input:
        mouse_ref = 'species/mouse/genome.fa',
        mouse_annotation='species/mouse/annotation.gtf',
        human_ref='species/human/genome.fa',
        human_annotation = 'species/human/annotation.gtf'
    output:
        human_index=directory('species/human/star_index'),
        mouse_index=directory('species/mouse/star_index')
    run:
        if not os.path.isdir('species/human/star_index'):
            shell("STAR --runMode genomeGenerate \
                 --runThreadN {threads} \
                 --genomeDir {output.human_index} \
                 --genomeFastaFiles {input.human_ref} \
                 --sjdbGTFfile {input.human_annotation}")

        if not os.path.isdir('species/mouse/star_index'):
            shell("STAR --runMode genomeGenerate \
                 --runThreadN {threads} \
                 --genomeDir {output.mouse_index} \
                 --genomeFastaFiles {input.mouse_ref} \
                 --sjdbGTFfile {input.mouse_annotation}")
        else:
            print('directory exists, continuing to alignment')

rule STAR_align:
    input:
        expand('{OUTDIR}/{sample}/mapping/{sample}_processed_{read}.fastq',OUTDIR=config['outdir'],sample=config['sample'],read=['R1','R2'])
    output:
        human=temp('{OUTDIR}/{sample}/mapping/human_Aligned_out.bam'),
        mouse=temp('{OUTDIR}/{sample}/mapping/mouse_Aligned_out.bam')
    params:
        human='species/human/star_index',
        mouse='species/mouse/star_index',
        human_annotation='species/human/annotation.gtf',
        mouse_annotation='speceies/mouse/annotation.gtf',
        aln="--readFilesCommand zcat --outSAMtype BAM Unsorted"
    threads:
        config['threads']
    shell:
        """
        STAR --genomeDir {params.human} \
            --readFilesIn {input} \
            --outFileNamePrefix {OUTDIR}/{sample}/mapping/human_ \
            --sjdbGTFfile {params.human_annotation} \
            {params.aln} \
            --runThreadN {threads}
        STAR --genomeDir {params.mouse} \
            --readFilesIn {input} \
            --outFileNamePrefix {OUTDIR}/{sample}/mapping/mouse_ \
            --sjdbGTFfile {params.mouse_annotation} \
            {params.aln} \
            --runThreadsN {threads}
        """


rule SortSam:
    input:
        human='{OUTDIR}/{sample}/mapping/human_Aligned_out.bam',
        mouse='{OUTDIR}/{sample}/mapping/mouse_Aligned_out.bam'
    output:
        human=temp('{OUTDIR}/{sample}/mapping/human_sorted_aligned.bam'),
        mouse=temp('{OUTDIR}/{sample}/mapping/mouse_sorted_aligned.bam')
    shell:
        """
        java -jar {picard} SortSam I={input.human} O={output.human} SORT_ORDER=queryname
        java -jar {picard} SortSam I={input.mouse} O={output.mouse} SORT_ORDER=queryname
        """

rule MergeAlignment:
    input:
        human_aln='{OUTDIR}/{sample}/mapping/human_sorted_aligned.bam',
        ubam='{OUTDIR}/{sample}/mapping/tagged_polyA_adapter_trimmed.bam',
        mouse_aln='{OUTDIR}/{sample}/mapping/mouse_sorted_aligned.bam'
    output:
        human=temp('{OUTDIR}/{sample}/mapping/human_merged.bam'),
        mouse=temp('{OUTDIR}/{sample}/mapping/mouse_merged.bam')
    params:
        human_ref='species/human/genome.fa',
        mouse_ref='species/mouse/genome.fa',
        other='INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false'
    shell:
        """
        java -jar {picard} MergeBamAlignment REFERENCE_SEQUENCE={params.human_ref}\
                                            UNMAPPED_BAM={input.ubam} \
                                            ALIGNED_BAM={input.human_aln} \
                                            OUTPUT={output.human} \
                                            {params.other}
        java -jar {picard} MergeBamAlignment REFERENCE_SEQUENCE={params.mouse_ref} \
                                            UNMAPPED_BAM={input.ubam} \
                                            ALIGNED_BAM={input.mouse_aln} \
                                            OUTPUT={output.mouse} \
                                            {params.other}
        """

rule TagReadWithGeneExon:
    input:
        human='{OUTDIR}/{sample}/mapping/human_merged.bam',
        mouse='{OUTDIR}/{sample}/mapping/mouse_merged.bam'
    output:
        human=temp('{OUTDIR}/{sample}/mapping/human_merged_tagged.bam'),
        mouse=temp('{OUTDIR}/{sample}/mapping/mouse_merged_tagged.bam')
    params:
        human_annotation='species/human/annotation.gtf',
        mouse_annotation='species/mouse/annotation.gtf'
    shell:
        """
        {dropseq}/TagReadWithGeneExon I={input.human} O={output.human} ANNOTATIONS_FILE={params.human_annotation} TAG=GE
        {dropseq}/TagReadWithGeneExon I={input.human} O={output.human} ANNOTATIONS_FILE={params.human_annotation} TAG=GE
        """

rule filer_mm_reads:
    input:
        human='{OUTDIR}/{sample}/mapping/human_merged.bam',
        mouse='{OUTDIR}/{sample}/mapping/mouse_merged_tagged.bam'
    output:
        human='{OUTDIR}/{sample}/mapping/human_final.bam',
        mouse='{OUTDIR}/{sample}/mapping/mouse_final.bam'
    shell:
        """
        python scripts/filter_mm_reads.py --in-bam {input.human} --out-bam {output.human}
        python scripts/filter_mm_reads.py --in-bam {input.mouse} --out-bam {output.mouse}
        """
