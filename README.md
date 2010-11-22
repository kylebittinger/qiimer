qiimer
======

This package provides R functions for (1) reading QIIME output files
and (2) creating figures from the resultant data frames.

Installation
------------

To install this library, type

    R CMD INSTALL /path/to/qiimer

where "/path/to/qiimer" is the path of the qiimer package directory
(the directory containing this file).

We have not set up a CRAN-style repository yet, so we are not yet able
to support install.packages() from the R shell.

Documentation
-------------

PDF documentation may be generated from embedded comments in the R
code.  The roxygen library is required for this process.  First,
generate the Rdoc files from the R source code by typing

    R CMD roxygen -d .

from the current directory.  Then, generate a PDF document by typing

    R CMD Rd2pdf .
