import os
import setuptools
from version import version as this_version


this_directory =  os.path.abspath(os.path.dirname(__file__))
version_path = os.path.join(this_directory, 'MoNeT_MGDrivE', '_version.py')
print(this_version)
with open(version_path, 'wt') as fversion:
    fversion.write('__version__ = "'+this_version+'"')

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='MoNeT_MGDrivE',
    version=this_version,
    url='https://github.com/Chipdelmal/MoNeT_MGDrivE',
    author='Hector M. Sanchez C.',
    author_email='sanchez.hmsc@berkeley.edu',
    description="MoNeT python package for MGDrivE data analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy', 'scipy', 'matplotlib', 'ipython',
        'jupyter', 'pandas', 'sympy', 'compress-pickle',
        'joblib', 'more-itertools'
    ],
    license='MIT',
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "Operating System :: OS Independent",
    ]
 )