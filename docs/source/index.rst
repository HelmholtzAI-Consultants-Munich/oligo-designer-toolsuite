.. oligo-designer-toolsuite documentation master file, created by
   sphinx-quickstart on Thu Oct 13 15:49:40 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


*Oligo Designer Toolsuite* - Lightweight Development of Custom Oligo Design Pipelines
=======================================================================================

Oligonucleotides (abbrev. oligos) are short, synthetic strands of DNA or RNA that are designed with respect to a specific target region and 
have many application areas, ranging from research to disease diagnosis or therapeutics. 
Oligos can be used as primers during DNA amplification, as probes for *in-situ* hybridization or as guide RNAs for CRISPR-based gene editing. 
Based on the intended application and experimental design, researchers have to customize the length, sequence composition, and thermodynamic 
properties of the designed oligos. 

|image0|

.. |image0| image:: _figures/oligo_design.png


Various tools exist that provide custom design of oligo sequences depending on the area of application. Interestingly, all those pipelines 
have many common basic processing steps, ranging from the generation of custom-length oligo sequences, the filtering of oligo sequences based on 
thermodynamic properties as well as the selection of an optimal set of oligos. Despite the fact that most tools apply the same basic processing 
steps, each newly developed tool usually uses its own implementation and different versions of package dependencies for those basic processing 
steps. As a consequence, the comparability of tools that differ only in certain steps is hampered, but also the maintenance of existing tools and the 
development of new tools is slowed down, because developers do not have a common resource for basic functionalities to use. We tackle 
this issue by providing such a common resource in our open-source *Oligo Designer Toolsuite*. 

**Oligo Designer Toolsuite is a collection of modules that provide all basic functionalities for custom oligo design pipelines within a flexible Python framework.** 
All modeles rely on a common underlying data structure, which allows the user to easily combine different modules, depending on the required processing steps.
We also provide ready-to-use oligo design pipelines for specific experimental setups, e.g. SCRINSHOT or SeqFISH+ probe design for Spatial Transcriptomics.



.. toctree::
   :maxdepth: 1
   :caption: GETTING STARTED

   _getting_started/installation.rst
   _getting_started/introduction_framework.rst
   _getting_started/run_ready_to_use_pipelines.rst


.. toctree::
   :maxdepth: 1
   :caption: TUTORIALS

   _tutorials/build_a_custom_pipeline

.. toctree::
   :maxdepth: 2
   :caption: API Documentation:

   oligo_designer_toolsuite.rst


Contributing
-------------

Contributions are more than welcome! Everything from code to notebooks to examples and documentation are all equally valuable so please don't feel you can't contribute. 
To contribute please fork the project make your changes and submit a pull request. We will do our best to work through any issues with you and get your code merged into the main branch.

How to cite
------------

If the Ologo Designer Toolsuite is useful for your research, consider citing the package:

..  code-block::

   @software{lisa_sousa_2023_7823048,
      author       = {  Isra Mekki, 
                        Francesco Campi, 
                        Louis Kümmerle,
                        Hanane Mohaouchane, 
                        Maksym Tretiakov, 
                        Anna Starovoit,
                        Cheng-Wei Liao,
                        Marie Piraud,
                        Lisa Barros de Andrade e Sousa},
      title        = {{Oligo Designer Toolsuite}},
      month        = april,
      year         = 2023,
      publisher    = {Zenodo},
      version      = {v0.1.3},
      doi          = {10.5281/zenodo.7823048},
      url          = {https://doi.org/10.5281/zenodo.7823048}
   }


License
--------

``oligo-designer-toolsuite`` is released under the MIT license. See `LICENSE <https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/blob/dev/LICENSE>`_ for additional details about it.
