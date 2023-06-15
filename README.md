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

# Usage

## Add Species Information
Required Inputs can be found at: https://www.ncbi.nlm.nih.gov/datasets/genome/
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

## Run Pipeline
```
snakemake --configfile workflow/config.yaml --cores <n_cores>
```
