from distutils.core import setup

setup(
    name='oligo-designer-toolsuite',
    version='0.1',
    packages=['oligo_designer_toolsuite'],
    install_requires=['datetime', 'argparse', 'pandas', 'iteration_utilities', 
                      'Bio', 'gtfparse', 'pyfaidx', 'pyyaml', 'pybedtools', 'networkx'],
    long_description=open('README.md').read(),
    author='Lisa Barros de Andrade e Sousa',
    entry_points={
        'console_scripts': [
            'padlock_probe_designer = oligo_designer_toolsuite.pipelines.padlock_probe_designer:main'
        ]
    },
)
