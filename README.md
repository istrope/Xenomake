# Xenomake
Processing pipeline for Patient Derived Xenograft (PDX) Spatial Transcriptomics Datasets

## PDX Model
Patient derived xenografts are a widely used component of cancer research <br>
These models consist of grafting a human tumor sample into mice to study interactions between a patient's tumor and stroma, screening compounds that target a patient's tumor, and simulating human tumor biology 

## Spatial Integration
With the rising popularity of spatially resolved transcriptomic technologies, there is a growing need for processing pipelines that can handle reads from PDX samples. Spatial processing has the potenntial to further analyze tumor/stroma interactions and tumor microenvironments. <br>

## Pipeline
To facilitate the adoption of spatial transcriptomics for PDX studies, we thus have developed a pipeline "Xenomake", which is an end-to-end pipeline that includes read alignment and gene quantification steps for reads generated by spatial transcriptomic platforms, and uses a xenograft sorting tool (Xengsort) to apportion xenograft reads to the host and graft genomes.<br>
Xenomake's structure is based on Snakemake and includes scripts to handle multi-mapping,downstream analysis, and handling "ambiguous" reads 


# Installation and Usage for Xeomake
check: https://github.com/Biivy/xenomake/wiki
