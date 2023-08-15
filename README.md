# Xenomake
Processing pipeline for Patient Derived Xenograft (PDX) Spatial Transcriptomics Datasets


# Build Xenomake

## Download Scripts
```
git clone https://github.com/Biivy/xenomake
cd xenomake
```

## Create conda environment and build dependencies
```
conda env create -n xenomake -f environment.yml
conda activate xenomake
```
## Install R Dependencies
```
install.packages(c('BiocManager','Matrix','data.table','R.utils'))
BiocManager::install(c('rhdf5','AnnotationDbi','org.Hs.eg.db','org.Mm.eg.db'))
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

# RUN Xeomake
check: https://github.com/Biivy/xenomake/wiki
