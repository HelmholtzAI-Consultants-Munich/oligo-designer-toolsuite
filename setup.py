from distutils.core import setup

setup(
    name='ProbeDesign',
    version='0.1',
    packages=['src'],
    install_requires=['argparse','pybedtools','gtfparse','Bio','pandas'],
    long_description=open('README.md').read(),
    author='Lisa B. A. Sousa & Erinc Merdivan',
)