# xenomake
Processing pipeline for Patient Derived Xenograft (PDX) Spatial Transcriptomics Datasets


# Build Xenomake

## Download Scripts
```
git clone https://github.com/Biivy/xenomake
cd xenomake
```

## Create conda environment and build dependencies
```
conda create -n xenomake -f environment.yaml
conda activate xenomake
```
## Install Xengsort for Xenograft Read Sorting <br>
see https://gitlab.com/genomeinformatics/xengsort for more information <br>
Dependencies for xengsort are included within the environment.yaml file
```
git clone https://gitlab.com/genomeinformatics/xengsort.git
cd xengsort  # the directory of the cloned repository
conda activate xenomake   # activate conda environment if not done already
pip install -e .
```
# Test Installation
#### Test Xengsort Installation
```
snakemake -s test/xengsort.smk --cores 4
```
#### Test Xenomake Pipeline
```
snakemake -s test_data/xenomake_test.smk --cores 4 -n  # Dry Run
snakemake -s test_data/xenomake_test.smk --dag  --cores 4| dot -Tsvg > dag.svg  #Directed acrylic graph
snakemake -s test_data/xenomake_test.smk  --cores 4 #Running through sample dataset
```
# Usage

## Add Species Information
Required Inputs can be found at: https://www.ncbi.nlm.nih.gov/datasets/genome/ <br>
Run <xenomake_dir>/scripts/species_parser.py --help for more information and all flag options

#### Required Flags:
1. **--mouse_ref:** Mouse Reference Assembly in fasta format
2. **--human_ref:** Human Reference Assembly in fasta format
3. **--mouse_annotation:** Mouse Genome Annotation File in gtf format
4. **--human_annotation:** Human Genome Annotation File in gtf format
#### Command Line Implementation:
```
python scripts/species_parser.py \
--mouse_ref <mouse_reference_assembly.fa> \
--human_ref <human_reference_assemble.fa> \
--mouse_annotation <mouse_annotation.gtf> \
--human_annotaion <human_annotation.gtf> \
```

## Create Configuration File
Use python <xenomake_dir>/scripts/config.py --help for more information and all flags
#### Required Flags:
1. **--r1:** paired end read in fastq format
2. **--r2:** second paired end read
3. **--outdir:** name of project directory where output file will be written
4. **--sample:** name of sample to prepend filenames
```
python scripts/config.py \
--r1 <sample_R1.fastq.gz> \
--r2 <sample_R2.fastq.gz> \
--outdir <dir_name> \
--sample <sample_name> \
--threads <n_threads>
```
#### Other Dependencies
Xenomake requires Dropseq-tools and picard to successfully run <br>
By default the path used by the pipeline is <cwd>/tools/Dropseq-tools-2.5.3 and <cwd>/tools/picard.jar <br>
However, you can specify the absolute path to your install with the **--picard** and **--dropseq_tools** flags in config.py <br>
Availability: https://github.com/broadinstitute/Drop-seq and https://broadinstitute.github.io/picard/
  
## Run Pipeline
```
snakemake -s xenomake.smk --configfile workflow/config.yaml --cores <n_cores>
```
