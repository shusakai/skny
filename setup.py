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

setup(
    author="Shunsuke A. Sakai",
    author_email='shusakai@east.ncc.go.jp',
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Spatial gene expression anlysis stratified by distance from tumor for multiple platform such as Xenium In Situ, CosMx, and PhenoCycler.",
    entry_points={
        'console_scripts': [
            'skny=skny.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='skny',
    name='skny',
    packages=find_packages(include=['skny', 'skny.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/shusakai/skny',
    version='0.1.1',
    zip_safe=False,
)
