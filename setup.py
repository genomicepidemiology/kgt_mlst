from setuptools import setup, find_packages

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

from kgt_mlst.version import __version__

setup(
    name='kgt_mlst',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version=__version__,
    packages=find_packages(),
    data_files=[],
    include_package_data=True,
    url='https://https://github.com/MBHallgren/kgt_mlst',
    license='',
    install_requires=(),
    author='Malte B. Hallgren',
    scripts=['bin/kgt_mlst'],
    author_email='malhal@food.dtu.dk',
    description='kgt_mlst',
)