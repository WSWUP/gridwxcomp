Installation
============

Currently we recommend using the provided conda environment file to install ``gridwxcomp`` and its dependencies in a virtual environment. Download the `environment.yml <https://raw.githubusercontent.com/WSWUP/gridwxcomp/master/gridwxcomp/env/environment.yml>`_ file and then install and activate it. If you don't have conda `get it here <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_. To install dependencies in a virtual environment run 

.. code-block:: bash

    $ conda env create -f environment.yml

To activate the environment before using ``gridwxcomp`` run

.. code-block:: bash

    $ conda activate gridwxcomp

After installing all the dependencies using conda, install ``gridwxcomp`` using `pip <https://pip.pypa.io/en/stable/installing/>`_,

.. code-block:: bash

    $ pip install gridwxcomp

Due to dependency conflicts you may have issues directly installing with pip before activating the conda environment. This is because the package includes several modules that are not pure Python such as GDAL and pyproj which seem to be better handled by conda. 

Alternatively, or if there are installation issues, you can manually install. First activate the ``gridwxcomp`` conda environment (above). Next, clone or download the package from `GitHub <https://github.com/WSWUP/gridwxcomp>`_ or `PyPI <https://pypi.org/project/gridwxcomp/>`_ and then install locally with pip in "editable" mode. For example with cloning,

.. code-block:: bash

    $ git clone https://github.com/WSWUP/gridwxcomp.git
    $ cd gridwxcomp

If you are experiencing errors on installing the ``gridwxcomp`` conda environment above with dependencies. For example, if the Shapely package is not installing from the enironment.yml file, remove it or modify it from the "setup.py" file in the install requirements section before you install gridwxcomp from source with:

.. code-block:: bash

    $ pip install -e .

More help with installation issues related to dependency conflicts can be found in the ``gridwxcomp`` `issues <https://github.com/WSWUP/gridwxcomp/issues>`_ on GitHub, be sure to check the closed issues as well.
