#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

#requirements = ['Click>=7.0', ]

test_requirements = [ ]

setup_requirements = []

setup(
    author="Shunsuke A. Sakai",
    author_email='shusakai@east.ncc.go.jp',
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Image Processing',
    ],
    description=(
        "SKNY (SpatialKNifeY): tumor microenvironment analysis for single-cell spatial omics "
        "(Xenium, CosMx, PhenoCycler). Tumor border detection, distance from the tumor, "
        "core / margin / distant regions, gene expression and cell density per region."
    ),
    long_description_content_type='text/x-rst',
    install_requires=requirements,
    extras_require={
        # stlearn is only needed to build the 10 um grid used by the expression-based workflow
        # (``stlearn.tl.cci.grid``); SKNY itself does not import it.
        'stlearn': ['stlearn'],
        'test': ['pytest'],
    },
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords=(
        "spatial transcriptomics, spatial omics, spatial proteomics, tumor microenvironment, "
        "tumor border, tumor margin, spatial domain, Xenium, CosMx, PhenoCycler, "
        "single-cell, AnnData, scanpy, squidpy, stlearn, cancer, bioinformatics"
    ),
    name='skny',
    packages=find_packages(include=['skny', 'skny.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/shusakai/skny',
    project_urls={
        'Documentation': 'https://skny.readthedocs.io',
        'Source': 'https://github.com/shusakai/skny',
        'Issue Tracker': 'https://github.com/shusakai/skny/issues',
        'Paper': 'https://doi.org/10.1371/journal.pcbi.1012854',
        'Changelog': 'https://github.com/shusakai/skny/blob/main/HISTORY.rst',
    },
    version='0.3.2',
    zip_safe=False,
)
