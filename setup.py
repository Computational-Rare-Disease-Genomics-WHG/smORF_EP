#!/usr/bin/env python

"""
Setup script for the Python package

# For development setup
pip install --editable .

# To build the package
python3 setup.py install

# Installing from Github 
pip3 install git+ssh://github.com/Computational-Rare-Disease-Genomics-WHG/smORF_EP.git

# To upload to pypi using twine 
python3 setup.py sdist bdist_wheel
twine upload dist/*
"""

import setuptools

with open("requirements.txt", encoding='utf-8') as f:
    requirements = f.read().splitlines()
    setuptools.setup(
        name='smorfep', # Don't change this name
        version='1.0.0',
        author='Computational Rare Disease Genomics WHG',
        description='smORF-EP tool',
        long_description=open('README.md', encoding='utf-8').read(),
        long_description_content_type='text/markdown',
        url='https://github.com/Computational-Rare-Disease-Genomics-WHG/smORF_EP',
        license='MIT',
        packages=['smorfep'],
        install_requires=requirements,
        include_package_data=True,
        zip_safe=False,
        # Add the installation installation script her
        scripts=['smorfep/download/compute_transcripts_gencode.py',  
                 'smorfep/download/preprocess_gff.py',
                 'smorfep/download/compute_introns_gencode_per_transc.py',
                 'smorfep/download/ref_per_chr.py', 
                 'smorfep/utils/genetic_code.py',
                 'smorfep/utils/functions.py',
                 'smorfep/utils/tool_script.py'
                 ], 
        entry_points={
            'console_scripts': [
                'smorfep = smorfep.run:main',
                'smorfinit = smorfep.downloader:main',
                'smorfinput = smorfep.input_generator:main',
                'checkseq = smorfep.check_ref_seq:main',
                'examplewizard = smorfep.generate_examples:main',
                'compare2vep = smorfep.compare2VEP:main'
            ],
        },
        keywords='bioinformatics',
        classifiers=[
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Unix',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    )