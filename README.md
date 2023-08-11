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
conda env create -n xenomake -f environment.yaml
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
Allow ~15 minutes to complete all tests <br>
If you just want to test default settings, run line two and three
```
snakemake -s test/xengsort.smk --cores 4  #all test parameters (440 jobs)
xengsort index --host test/xengsort_data/ecoli.fa.gz --graft test/xengsort_data/ehec.fa.gz -n 10000 --index test/test_index.zarr
xengsort classify --index test/test_index.zarr --fastq test/xengsort_data/gzipped_chunks1.fq.gz --out test/test_classify
```
#### Test Xenomake Pipeline
Allow ~1hr to complete <br> <br>
#Download reference genome assemblies
```
#Download mm10 genome
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001635.27/download?include_annotation_type=GENOME_FASTA,GENOME_GTF&filename=GCF_000001635.27.zip" -H "Accept: application/zip"
mv ncbi_dataset mm10/
#Download hg38 genome
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA,GENOME_GTF&filename=GCF_000001405.40.zip" -H "Accept: application/zip"
mv ncbi_dataset hg38
```
```
chmod +x test/xenomake_test.sh
bash test/xenomake_test.sh
```
