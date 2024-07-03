.. _faq:

FAQ
===

How are monthly station to gridded data bias ratios calculated?
   The following steps describes the default method to the ``calc_bias_ratios`` function (i.e., the ``method`` kwarg = 'long_term_mean') for all variables except temperature:

   #. Pair up daily station data with corresponding gridded data.
   #. For each month on record, select those that have at least 10 (default threshold) days of paired data. 
   #. Group all daily data by month, e.g. all dates that lie in the month of July.
   #. Take the mean of the grouped station data and divide by the mean of the grouped gridMET data.

       * For temperature data the difference is calculated as opposed to station:gridMET ratios. In other words step 4. becomes: mean of grouped station minus the mean of grouped gridMET data.

