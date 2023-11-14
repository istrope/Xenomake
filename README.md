1. [About Xenomnake](#sec1) </br>
    1.1. [PDX Model](#sec1.1)</br>
    1.2. [Spatial Integration](#sec1.2)</br>
    1.3. [Pipeline Overview](#sec1.3)</br>
    1.4. [Publication](#sec1.4)</br>
2. [Installation](#sec2)</br>
    2.1. [Download Xenomake Scripts](#sec2.1)</br>
    2.2. [Create conda environment and build dependencies](#sec2.2)</br>
    2.3. [Install Xengsort for xenograft read sorting](#sec2.3)</br>
    2.4. [Install genomes and indices](#sec2.4)</br>
    2.5. [Download optional tools](#sec2.4)</br>
3. [Running Xenomake](#sec3)</br>
	3.1. [QuickStart](#sec3.1)<br>
 	3.1.2 [Memory Requirements and Runtimes](#sec3.1.2)
    3.2. [Step1: Add Species Information](#sec3.2)</br>
    3.2.1. [Description](#sec3.2.1)<br>
		3.2.2. [Flags](#sec3.2.2)<br>
		3.2.3. [Command Line Implementation](#sec3.2.3)<br>
		3.2.4. [Output](#sec3.2.4)<br>
    3.3. [Step2: Project Configuration](#sec3.3)</br>
    	3.3.1 [Description](#sec3.3.1)<br>
		3.3.2 [Flags](#sec3.3.2)<br>
		3.3.3. [Command Line Implementation](#sec3.3.3)<br>
		3.3.4. [Output](#sec3.3.4)<br>
    3.4. [Step3: Run Snakemake Pipeline](#sec3.4)</br>
		3.4.1. [Snakemake Description](#sec3.4.1)<br>
		3.4.2. [Dry Run](#sec3.4.2)<br>
		3.4.3. [Execute Snakemake Workflow](#sec3.4.3)<br>
		3.4.4. [Snakemake Standard Outs](#sec3.4.4)<br>
		3.4.5. [Xenomake Output Structure](#sec3.4.4)<br>
4. [Downstream Cell-Cell Interactions](#sec4)<br>
5. [QC](#sec5)</br>
    5.1. [Scanpy QC Plots](#sec5.1)</br>
    5.2. [Xengsort Metrics](#sec5.2)</br>
    5.3. [STAR Metrics](#sec5.3)</br>
    5.4. [UMI Metrics](#sec5.4)</br>
6. [Test Dataset](#sec6)</br>
	6.1. [Download Data](#sec6.1)</br>
	6.2. [Run Xenomake Pipeline](#sec6.2)</br>
       
<a name="sec1"></a>
# About Xenomake
Processing pipeline for Patient Derived Xenograft (PDX) Spatial Transcriptomics Datasets

<a name="sec1.1"></a>
## PDX Model
Patient derived xenografts are a widely used component of cancer research <br>
These models consist of grafting a human tumor sample into mice to study interactions between a patient's tumor and stroma, screening compounds that target a patient's tumor, and simulating human tumor biology 
<a name="sec1.2"></a>
## Spatial Integration
With the rising popularity of spatially resolved transcriptomic technologies, there is a growing need for processing pipelines that can handle reads from PDX samples. Spatial processing has the potenntial to further analyze tumor/stroma interactions and tumor microenvironments. <br>
<a name="sec1.3"></a>
## Pipeline
To facilitate the adoption of spatial transcriptomics for PDX studies, we thus have developed a pipeline "Xenomake", which is an end-to-end pipeline that includes read alignment and gene quantification steps for reads generated by spatial transcriptomic platforms, and uses a xenograft sorting tool (Xengsort) to apportion xenograft reads to the host and graft genomes.<br>
Xenomake's structure is based on Snakemake and includes scripts to handle multi-mapping,downstream analysis, and handling "ambiguous" reads 
<a name="sec1.4"></a>
## Publication
Xenomake: a pipeline for processing and sorting xenograft reads from spatial transcriptomic experiments <br>
Benjamin S Strope, Katherine E Pendleton, William Z Bowie, Gloria V Echeverria, Qian Zhu <br>
bioRxiv 2023.09.04.556109; doi: https://doi.org/10.1101/2023.09.04.556109
<a name="sec2"></a>
# Installation
<a name="sec2.1"></a>
## Download Scripts
```
git clone https://github.com/Biivy/Xenomake
# Move into repository
cd Xenomake
chmod -R +x scripts/
```
<a name="sec2.2"></a>
### Create conda environment and build dependencies
```
# Create conda environment with specified requirements
conda env create -n xenomake -f environment.yaml
# Activate conda environment
conda activate xenomake
```
<a name="sec2.3"></a>
### Install Xengsort for Xenograft Read Sorting
see https://gitlab.com/genomeinformatics/xengsort for more information
Dependencies for xengsort are included within the environment.yaml file
```
git clone https://gitlab.com/genomeinformatics/xengsort.git
cd xengsort  # the directory of the cloned repository
conda activate xenomake   # activate conda environment if not done already
pip install -e .
```
<a name="sec2.4"></a>
## Install Genomes and Indices

### Zenodo Repository
Contains Genome Files, Annotation Files, and Xengsort Index Files <br>
Mouse Version: mm10 <br>
Human Version: hg38 <br>
Zenodo Repository: https://zenodo.org/record/8250705

### NCBI Genomes 
Optional in case users want to use a different assembly version <br>
Webpage: https://www.ncbi.nlm.nih.gov/datasets/genome/ 

<a name="sec2.5"></a>
## Download Dropseq (optional)
Current version of Dropseq-tools is 2.5.3 and is present at <repo_dir>/tools/Dropseq-2.5.3/ <br>
For updated versions check:
https://github.com/broadinstitute/Drop-seq 
## Download Picard (optional)
### Version Control
Currently, the provided picard.jar file is compatible with java installed in Xenomake conda environment. <br>
Java: v11.0.15 <br>
Picard: v2.27.5 <br>

### Picard Github
If you choose to use your own picard executable, ensure compatibility with java <br>
https://broadinstitute.github.io/picard/
<a name="sec3"></a>
# Running Xenomake

<a name="sec3.1"></a>
## Quickstart
### Quickstart: Add Species Information:
```
# Print Help Function Call
python <repo_dir>/scripts/species_parser.py --help

# Create species folder 
python <reop_dir>/scripts/species_parser.py \
--mouse_ref <mouse_reference_assembly.fa> \
--human_ref <human_reference_assemble.fa> \
--mouse_annotation <mouse_annotation.gtf> \
--human_annotaion <human_annotation.gtf> \
```


### Quickstart: Create Configuration File:
```
# Print Help Function Call
python <repo_dir>/scripts/config.py --help

# Create config.yaml 
python <reop_dir>/scripts/config.py \
--repository <path_to_Xenomake_dir>
--r1 <sample_R1.fastq.gz> \
--r2 <sample_R2.fastq.gz> \
--outdir <dir_name> \
--sample <sample_name> \
--threads <n_threads>
```


### Quickstart: Run Snakemake Pipeline
```
# Dry Run
snakemake -s <repo_dir>/snakemake/main.smk -n --cores 8
# Initialize Xenomake Pipeline (species dir and config.yaml need to be created)
snakemake -s <repo_dir>/snakemake/main.smk --cores <n_threads> --keep-going
```
<a name="sec3.1.2"></a>
### Memory Requirements and Runtimes
|               Process               | Threads |        Time        | RAM Recommended |
|:-----------------------------------:|:-------:|:------------------:|:---------------:|
|            Xengsort Index           |    8    |       33 min       |      25 Gb      |
|              STAR Index             |    8    | ~90 min per genome |      128 Gb     |
| Xenomake Pipeline<br>(Test Dataset) |    8    |       37 min       |      35 Gb      |
| Xenomake Pipeline<br>(Test Dataset) |    4    |       42 min       |      35 Gb      |

<a name="sec3.2"></a>
## Step 1: Add Species Information:
This is a **REQUIRED** step to initialize your implementation of the xenomake pipeline
<a name="sec3.2.1"></a>
### Description
species_parser.py links all inputs in a new directory structure that is then easily referenced by snakemake
<a name="sec3.2.2"></a>
### Required Flags: 
- **--mouse_ref:** Mouse Reference Assembly in fasta format<br>
- **--human_ref:** Human Reference Assembly in fasta format<br>
- **--mouse_annotation:** Mouse Genome Annotation File in gtf format<br>
- **--human_annotation:** Human Genome Annotation File in gtf format<br>
### Optional Flags:
Including optional flags for indices will skip parts in the pipeline such as *STAR --runMode CreateIndex* and *xengsort index*<br>
These two rules can be time costly and if these files are already generated, it doesn't need to be done again. <br> <br>
- **--human_index:** STAR index directory for human genome reference<br>
- **--mouse_index:** STAR index directory for mouse genome reference<br>
- **--xengsort_index:** Xengsort Index output to be used in xengsort classify step
<a name="sec3.2.3"></a>
### Command Line Implementation:
```
# Print Help Function Call
python scripts/species_parser.py --help

# Create species folder 
python <reop_dir>/scripts/species_parser.py \
--mouse_ref <mouse_reference_assembly.fa> \
--human_ref <human_reference_assemble.fa> \
--mouse_annotation <mouse_annotation.gtf> \
--human_annotaion <human_annotation.gtf> \
```
<a name="sec3.2.4"></a>
### Output Message
```
species directory successfully completed, keep project execution in the current directory
```

### Output Directory Structure
**├─ species** <br>
│   ├── idx.hash <br>
│   ├── idx.info <br>
**│   ├── human** <br>
│   │   ├── genome.fa <br>
│   │   ├── annotation.gtf <br>
│   │   ├── star_index/ <br>
**│   ├── mouse** <br>
│   │   ├── genome.fa <br>
│   │   ├── annotation.gtf <br>
│   │   ├── star_index/ 
<a name="sec3.3"></a>
## Step 2: Project Configuration
This is a **REQUIRED** step to initialize your implementation of the xenomake pipeline
<a name="sec3.3.1"></a>
### Description
config.py generates a .yaml file that is referenced by snakemake that defines locations of inputs, dependencies, filenames, and output directory
<a name="sec3.3.2"></a>
### Required Flags:
- **--repository**: path to cloned Xenomake Directory
- **--r1:** paired end read in fastq format <br>
- **--r2:** second paired end read <br>
- **--outdir:** name of project directory where output file will be written <br>
- **--sample:** name of sample to prepend filenames
- **--repository**:path to the xenomake cloned repository on home system

### Optional Flags:
- **--threads**: number of cores to use<br>
- **--picard**: path to user preferred picard.jar tool<br>
- **--dropseq_tools**: path to user preferred download of dropseq tools directory<br>
- **--downstream**: performs downstream processing using scanpy (default: True) <br>
- **--genic_only**: only uses UTR and Exonic reads for umi matrix (default: True) <br>
- **--mouse_assembly**: version of mouse genome used (default: mm10) <br>
- **--human_assembly**: version of human genome used (default: hg38) <br>
- **--barcode**: list of visium barcodes (default: visium v2) <br>
- **--ambiguous**: handle multimap from "both" and "ambiguous" output reads (default: True)
<a name="sec3.3.3"></a>
### Command Line Implementation:
```
# Print Help Function Call
python scripts/config.py --help

# Create config.yaml 
python <reop_dir>/scripts/config.py \
--repository <path_to_Xenomake_dir>
--r1 <sample_R1.fastq.gz> \
--r2 <sample_R2.fastq.gz> \
--outdir <dir_name> \
--sample <sample_name> \
--threads <n_threads>
```
<a name="sec3.3.4"></a>
### Output Message
Prints When config.py completes successfully
```
assembly versions ['hg38', 'mm10']
downstream processing: True
ambiguous reads handling: True
cores: 8
project "downsampled" initialized, proceed to snakemake execution
```

### Example Configuration File
![config.yaml](https://github.com/istrope/Xenomake/blob/main/figures/config.png)
<a name="sec3.4"></a>
## Step 3: Run Snakemake Pipeline
<a name="sec3.4.1"></a>
### Snakemake
Xenomake uses the Snakemake workflow management system to create reproducible and scalable data analyses <br>
To get an understanding of snakemake, please visit: https://snakemake.github.io/ <br>
For additional information, read Snakemake's paper https://doi.org/10.12688/f1000research.29032.1

<a name="sec3.4.2"></a>
### Dry Run
It is recommended to first execute a dry-run with flag **-n**, Snakemake will only show the execution plan instead of performing steps. <br>
### Requirements: <br>
**-s**: snakemake file containing final target files <br>
by default, snakemake searches for a Snakemake.smk file, but for Xenomake the target files are contained within main.smk <br>
**-n**: Dry-Run showing the workflow <br>
**--cores**: Number of cores to run Xenomake workflow with (**required**)
```
# Dry Run
snakemake -s <source_dir>/snakemake/main.smk -n --cores 8
```
<a name="sec3.4.3"></a>
### Execute Xenomake Workflow
note: species directory and config.yaml need to be created from previous steps <br>
### Flags
**-s**: Snakemake target file <br>
**--cores**: number of cores used for execution <br>
**--keep-going**: Snakemake flag that will continue the workflow when possible if a single job fails but is not required for subsequent steps
```
snakemake -s <source_dir>/snakemake/main.smk --cores <n_threads> --keep-going
```

<a name="sec3.4.4"></a>
### Example Snakemake Standard Outs
Below shows an example of the messages printed by Snakemake. Briefly, they give an overview of the jobs in line for execution, and print messages when a rule is being performed. This includes a few important features: <br> <br>

1. Start time <br>
2. Finish time <br>
3. Rule name <br>
4. Input files <br>
5. Output file destinations <br>
6. Wildcard values <br>
7. Additional params and log destinations <br> <br>

Snakemake will also print any standard out messages from tools performed. We have saved some of these as log files for subsequent review

```
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 16
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	DGE_Human
	1	Expression_QC
...
	1	Subset_Overlapped
	1	Xengsort_Clasify
	1	all
	28
Select jobs to execute...

[Fri Sep  8 10:36:04 2023]
rule Preprocess:
    input: fastq/sub1.fq, fastq/sub2.fq
    output: out/treated/preprocess/unaligned_bc_umi_tagged.bam, out/treated/logs/Preprocess.summary
    log: out/treated/logs/Preprocess.log
    jobid: 6
    wildcards: OUTDIR=out, sample=treated

[Fri Sep  8 10:36:06 2023]
Finished job 6.
1 of 28 steps (4%) done
Select jobs to execute...

...
...
...

[Fri Sep 8 10:47:41 2023]
localrule all:
    input: out/treated/final/treated_human_counts.tsv.gz, out/treated/final/treated_mouse_counts.tsv.gz, out/treated/final/treated_human.hdf5, out/treated/final/treated_mouse.hdf5, out/treated/final/treated_human.h5ad, out/treated/final/treated_mouse.hdf5, out/treated/final/treated_mouse.h5ad, out/treated/final/treated_human_processed.h5ad, out/treated/final/treated_mouse_processed.h5ad, out/treated/qc/human_qc.png
    jobid: 0
    threads: 16

[Fri Sep 8 14:37:01 2023]
Finished job 0.
29 of 29 steps (100%) done

```

<a name="sec3.4.5"></a>
### Output Directory Structure
├─ **OUTDIR** <br>
**│   ├── sample**<br>
│   │   ├── final.bam:  Mapped, Sorted, Xenograft Processed alignments in bam format <br>
**│   │   ├── final**<br>
│   │   │   ├── tissue_positions_list.csv<br>
│   │   │   ├── {sample}_counts.tsv.gz: *Gene by Barcode umi count matrix generated by dropseq DigitalGeneExpression*  <br>
│   │   │   ├── {sample}_human.h5ad:*Spatial Anndata object*  <br>
│   │   │   ├── {sample}_human.hdf5: *10X genomics formatted output structure (can be used for downstream spatial processing)*  <br>
**│   ├── xengsort:** *Xenograft sorted fastq files generated using xengsort toolkit* <br>
│   │   │   ├── ambiguous.fq.gz <br>
│   │   │   ├── both.fq.gz <br>
│   │   │   ├── graft.fq.gz <br>
│   │   │   ├── host.fq.gz <br>
│   │   │   ├── neither.fq.gz <br>
**│   ├── mapping** <br>
│   │   │   ├── STAR aligned genomes: *Aligned files without Xenograft Processing* <br>
**│   ├── logs** <br>
│   │   │   ├── *Tagging and trimming logs*<br>
│   │   │   ├── *STAR log files*<br>
│   │   │   ├── *BAM subsetting logs*<br>
│   │   │   ├── *Xengsort log*<br>
│   │   │   ├── *Multi-map filtering logs*<br>
│   │   │   ├── *Digital Gene Expression logs*<br>
**│   ├── qc** <br>
│   │   │   ├── *Xengsort summary files*<br>
│   │   │   ├── *STAR summary files*<br>
│   │   │   ├── *Dataset qc plots*<br>
│   │   │   ├── *Digital Gene Expression summary*<br>
**│   ├── preprocessing** <br>
│   │   │   ├── unaligned_bc_umi_tagged.bam: Unaligned bam file. Cell-barcode/UMI tagged and Adapter/PolyA trimmed using dropseq workflows <br>
<a name="sec4"></a>
# Downstream Cell-Cell Interactions
## Intro
Using the outputs from Xenomake, we were able to identify group specific expression and differential signalling genes in our model. We provide scripts in the repository `CellInteract` to perform identical analyses that are present in the Xenomake Publication
## Cell-Cell Signalling
Cell Signalling molecules such as chemokines/cytokines are simultaneously secreted by both stromal cell types and epithelial tumor cells in the tumor microenvironment. It is often difficult to attribute signalling expression to stromal or tumor cells in standard models. However in PDX models processed with means similar to Xenomake, we are able to study the differential production of cellular signals and identify biomarkers in stromal and epithelial compartments
## Xenomake Publication Figure (Organism Specific Signalling)
![Xenomake Plot](https://github.com/istrope/Xenomake/blob/main/figures/human_mouse_split.jpg)<br>
## Repository
Access further description as well as scripts to perform these downstream analysis:
`Link to CelllInteract Repository`: https://github.com/bernard2012/CellInteract
<a name="sec5"></a>
# QC
<a name="sec5.1"></a>
## Scanpy QC Plots
![QC Metrics](https://github.com/istrope/Xenomake/blob/main/figures/scanpy_qc.png) <br>
<a name="sec5.2"></a>
## Xengsort Metrics
Counts number of reads partitioned to each file

| prefix | host    | graft    | ambiguous | both    | neither |
|--------|---------|----------|-----------|---------|---------|
| test   | 6652723 | 13654867 | 2603727   | 1297366 | 154218  |
```
## Running time statistics
- Running time in seconds: 225.0 sec
- Running time in minutes: 3.75 min
- Done. [2023-08-17 16:19:49]
```

<a name="sec5.3"></a>
## STAR Metrics
```
                                 Started job on |       Aug 16 18:57:26
                             Started mapping on |       Aug 16 18:58:54
                                    Finished on |       Aug 16 19:04:22
       Mapping speed, Million of reads per hour |       968.48

                          Number of input reads |       88239353
                      Average input read length |       108
                                    UNIQUE READS:
                   Uniquely mapped reads number |       36755354
                        Uniquely mapped reads % |       41.65%
                          Average mapped length |       109.05
                       Number of splices: Total |       7555478
            Number of splices: Annotated (sjdb) |       7421066
                       Number of splices: GT/AG |       7487922
                       Number of splices: GC/AG |       23361
                       Number of splices: AT/AC |       1030
               Number of splices: Non-canonical |       43165
                      Mismatch rate per base, % |       1.75%
                         Deletion rate per base |       0.04%
                        Deletion average length |       1.90
                        Insertion rate per base |       0.04%
                       Insertion average length |       1.89
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       10085928
             % of reads mapped to multiple loci |       11.43%
        Number of reads mapped to too many loci |       768576
             % of reads mapped to too many loci |       0.87%
                                  UNMAPPED READS:
  Number of reads unmapped: too many mismatches |       0
       % of reads unmapped: too many mismatches |       0.00%
            Number of reads unmapped: too short |       38624549
                 % of reads unmapped: too short |       43.77%
                Number of reads unmapped: other |       2004946
                     % of reads unmapped: other |       2.27%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%
```
<a name="sec5.4"></a>
## UMI Metrics
Example of top 4 spots with highest expression
| CELL_BARCODE     | NUM_GENIC_READS | NUM_TRANSCRIPTS | NUM_GENES |
|------------------|-----------------|-----------------|-----------|
| CAGTTCCGCGGGTCGA | 429279          | 67931           | 9287      |
| GAAACTCGTGCGATGC | 398339          | 59154           | 8748      |
| AAGAGATGAATCGGTA | 362017          | 57447           | 8723      |
| CCAAGAAAGTGGGCGA | 367515          | 54742           | 8306      |
<a name="sec6"></a>
# Run Xenomake on Medulloblastoma Test Dataset
Follow Installation Instructions for Xenomake if not done already <br>
Perform all subsequent steps within Xenomake Directory
<a name="sec6.1"></a>
## Dataset
Paired End Fastq Files Downsamples to 8.4 million reads <br>
Taken from: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01185-4 <br>
### Download Dataset
```
wget https://zenodo.org/records/10120014/files/sub1.fq.gz
wget https://zenodo.org/records/10120014/files/sub2.fq.gz
```
### Download Genomes
```
wget https://zenodo.org/record/8326731/files/hg38_assembly.zip
unzip hg38_assembly.zip
wget https://zenodo.org/record/8326731/files/mm10_assembly.zip
unzip mm10_assembly.zip
```
<a name="sec6.2"></a>
## Run Xenomake Pipeline
### Activate Xenomake Environment
```
conda activate xenomake
```
### Species Setup
```
python scripts/species_parser.py \
--human_ref ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna \
--human_annotation ncbi_dataset/data/GCF_000001405.40/genomic.gtf \
--mouse_ref ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna \
--mouse_annotation ncbi_dataset/data/GCF_000001635.27/genomic.gtf
```
### Configuration
```
python scripts/config.py \
--repository ./ \
--r1 test/sub1.fq \
--r2 test/sub2.fq \
--outdir out \
--sample test
```
### Run Pipeline
Allow time for Xengsort Indexing and STAR indexing (This is the longest part)
```
snakemake -s snakemake/main.smk --cores 8
```
