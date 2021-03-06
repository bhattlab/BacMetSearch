from setuptools import setup, find_packages

setup(
    name="bacmetsearch",
    version='0.0.1_dev',
    description='A command line tool to perform local DIAMOND searches against the BacMet database.',
    url='https://github.com/bhattlab/BacMetSearch',
    author="Matt Durrant",
    author_email="mdurrant@stanford.edu",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'wget==3.2',
        'click==7.0',
        'biopython==1.76'
    ],
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'bacmetsearch = bacmetsearch.main:cli',
        ],
}
)
