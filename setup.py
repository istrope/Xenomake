from setuptools import setup

setup(
    #metadata
    name='xenomake',
    version='1.0.0',
    description='Xenomake: Processing Pipeline for Patient Derived Xenograft (PDX) Spatial Transcriptomics Datasets',
    url='https://github.com/istrope/Xenomake/',
    author='Ivy Strope',
    author_email='ivystrope@gmail.com',
    packages=['xenomake'],
    install_requires=['snakemake','argparse','pyyaml','scanpy'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3'
    ],
    
    #options
    python_requires='>=3.10',
    package_dir={
        'xenomake':'xenomake'
    },

    #package data
    package_data = {
        'xenomake': ['snakemake/*.smk',
                     'scripts/*.py',
                     'scripts/*.sh',
                     'data/barcodes/*.txt',
                     'data/barcodes/*.csv',
                     'data/Drop-seq_tools-2.5.3/*',
                     'data/urls/*.txt',
                     'data/picard.jar',
                     'data/*.yaml']
    },

    #entry points
    entry_points = {
        'console_scripts': [
            'xenomake = xenomake.smk:cmdline',
        ],
    },
)
