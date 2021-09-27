from distutils.core import setup

setup(
    name='ProbeDesign',
    version='0.1',
    packages=['src'],
    install_requires=['datetime', 'logging', 'os', 'multiprocessing', 'itertools', 'time', 'argparse', 'yaml', 'gzip', 'shutil', 'ftplib', 'pandas', 'pybedtools', 'gtfparse', 'pyfaidx', 'Bio'],
    long_description=open('README.md').read(),
    author='Lisa Barros de Andrade e Sousa & Erinc Merdivan',
)