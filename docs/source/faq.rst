FAQ
===

How are monthly station to gridded data bias ratios calculated?
   The following steps describes the default method to the :func:`gridwxcomp.calc_bias_ratios.calc_bias_ratios`
   function (i.e., the ``method`` kwarg = 'long_term_mean') for all variables
   except temperature:

   #. Pair up daily station data with corresponding gridded data.
   #. For each month on record, select those that have at least 10 (default threshold) days of paired data. 
   #. Group all daily data by month, e.g. all dates that lie in the month of July.
   #. Take the mean of the grouped station data and divide by the mean of the grouped gridded data.

       * For temperature data the difference is calculated as opposed to station:gridded ratios. In other words step 4. becomes: mean of grouped station minus the mean of grouped gridded data.

.. seealso::
    For more information on this calculation and the other options see the
    :ref:`Step 3: Calculate monthly, seasonal, and annual station:gridded
    biases and statistics` section in the tutorial. 
