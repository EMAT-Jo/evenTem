from setuptools import setup, find_packages

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md')) as f:
    long_description = f.read()


setup(
    name='eventem',
    version='0.0.1',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib'
    ],
    include_package_data=True,
    package_data={
        '': ['*/*.pyd','*/*.dll','*/*.so'],
    },
)