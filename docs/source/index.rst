.. oligo-designer-toolsuite documentation master file, created by
   sphinx-quickstart on Thu Oct 13 15:49:40 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*Oligo Designer Toolsuite* - Lightweight Development of Custom Oligo Design Pipelines
=======================================================================================


.. include:: ../../README.md
   :parser: myst_parser.sphinx_
   :start-after: <!-- LINK INTRODUCTION START -->
   :end-before: <!-- LINK INTRODUCTION END -->

.. toctree::
   :maxdepth: 1
   :caption: GETTING STARTED

   _getting_started/installation.rst
   _getting_started/introduction_framework.rst
   _getting_started/build_a_custom_pipeline.rst

.. toctree::
   :maxdepth: 1
   :caption: PIPELINES

   _pipelines/genomic_region_generator
   _pipelines/scrinshot_probe_designer
   _pipelines/merfish_probe_designer
   _pipelines/seqfishplus_probe_designer
   _pipelines/oligoseq_probe_designer


.. toctree::
   :maxdepth: 2
   :caption: API Documentation:

   _api_docs/oligo_designer_toolsuite.rst


Contributing
-------------

Contributions are more than welcome! Everything from code to notebooks to examples and documentation are all equally valuable so please don't feel you 
can't contribute. To contribute please fork the project make your changes and submit a pull request. We will do our best to work through any issues with 
you and get your code merged into the main branch.

For any further inquiries please send an email to `Lisa Barros de Andrade e Sousa <mailto:lisa.barros@helmholtz-munich.de>`_
or `Isra Mekki <mailto:isra.mekki@helmholtz-munich.de>`_.

How to cite
------------

If the Oligo Designer Toolsuite is useful for your research, consider citing the package:

::

   @software{campi_2023_7823048,
      author       = { Isra Mekki,
                     Francesco Campi,  
                     Louis Kümmerle,
                     Chelsea Bright,
                     Malte Lücken
                     Fabian Theis,
                     Marie Piraud,
                     Lisa Barros de Andrade e Sousa
                     },
      title        = {{Oligo Designer Toolsuite}},
      year         = 2023,
      publisher    = {Zenodo},
      version      = {v0.1.3},
      doi          = {10.5281/zenodo.7823048},
      url          = {https://doi.org/10.5281/zenodo.7823048}
   }


License
--------

``oligo-designer-toolsuite`` is released under the MIT license. See `LICENSE <https://github.com/HelmholtzAI-Consultants-Munich/oligo-designer-toolsuite/blob/dev/LICENSE>`_ for additional details about it.
