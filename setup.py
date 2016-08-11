"""

See:
https://github.com/Illumina/tHapMix
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path, walk

def gen_data_files(*dirs):
    results = []
    for src_dir in dirs:
        for root, dirs, files in walk(src_dir):
            results.append((root, map(lambda f:root + "/" + f, files)))
    return results
    
here = path.abspath(path.dirname(__file__))

setup(
    name='tHapMix',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='1.0',

    description='tHapMix',
    long_description='tHapMix - lightweight and fast simulation of tumour whole-genome sequencing datat',

    # The project's main homepage.
    url='https://github.com/Illumina/tHapMix',

    # Author details
    author='Sergii Ivakhno',
    author_email='pypa-dev@googlegroups.com',

    # Choose your license
    license='GPLv3',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 5 - Stable',

        # Pick your license as you wish (should match "license" above)
        'License :: GPLv3',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
  
    install_requires=['pysam', 'bx-python', 'HTSeq'],
    
    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=['bin'],
    
    data_files = gen_data_files("example", "config", "redist", "scripts", "pyflow"),
    
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'sample=sample:main',
        ],
    },
)