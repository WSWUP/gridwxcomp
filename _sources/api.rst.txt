API Reference
=============

This page documents all routines provided by ``gridwxcomp``. Tools for different routines are listed roughly in the order that they would normally be used.


.. contents:: :local:


Command Line Interface
^^^^^^^^^^^^^^^^^^^^^^

.. click:: gridwxcomp.scripts.gridwxcomp:gridwxcomp
  :prog: gridwxcomp

.. click:: gridwxcomp.scripts.gridwxcomp:prep_input
  :prog: gridwxcomp prep-input

.. click:: gridwxcomp.scripts.gridwxcomp:calc_bias_ratios
  :prog: gridwxcomp calc-bias-ratios

.. click:: gridwxcomp.scripts.gridwxcomp:download_gridmet_opendap
  :prog: gridwxcomp download-gridmet-opendap

.. click:: gridwxcomp.scripts.gridwxcomp:spatial
  :prog: gridwxcomp spatial

.. click:: gridwxcomp.scripts.gridwxcomp:plot
  :prog: gridwxcomp plot

Python functions, classes, and modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

prep\_input
-----------

.. autofunction:: gridwxcomp.prep_input


download\_gridmet\_opendap
--------------------------

.. autofunction:: gridwxcomp.download_gridmet_opendap


calc\_bias\_ratios
------------------

.. automodule:: gridwxcomp.calc_bias_ratios
    :members:
    :exclude-members: main, arg_parse
    :undoc-members:
    :show-inheritance:

spatial
-------

.. automodule:: gridwxcomp.spatial
    :members:
    :exclude-members: main, arg_parse, get_cell_ID
    :undoc-members:
    :show-inheritance:


plot
----

.. automodule:: gridwxcomp.plot
    :members:
    :exclude-members: arg_parse
    :undoc-members:
    :show-inheritance:

InterpGdal
----------

.. autoclass:: gridwxcomp.InterpGdal
    :members:
    :undoc-members:
    :show-inheritance:

util
----

.. automodule:: gridwxcomp.util
    :members:
    :undoc-members:
    :show-inheritance:

