#!/usr/bin/env python

"""
Setup script for the Python package
- Used for development setup with `pip install --editable .`
- Parsed by conda-build to extract version and metainfo
"""

import setuptools

PKG = 'smORF-EP'

with open("requirements.txt", encoding='utf-8') as f:
    requirements = f.read().splitlines()
    setuptools.setup(
        name=PKG,
        version='1.0.0',
        author='Computational Rare Disease Genomics WHG',
        description='My python project',
        long_description=open('README.md', encoding='utf-8').read(),
        long_description_content_type='text/markdown',
        url='https://github.com/Computational-Rare-Disease-Genomics-WHG/smORF_EP',
        license='MIT',
        packages=['smorfep'],
        install_requires=requirements,
        include_package_data=True,
        zip_safe=False,
        scripts=['smORF_EP.py'],
        keywords='bioinformatics',
        classifiers=[
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX',
            'Operating System :: Unix',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    )