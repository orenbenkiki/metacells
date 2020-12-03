Metacells - Single-cell RNA Sequencing Analysis
===============================================

.. image:: https://travis-ci.org/tanaylab/metacells.svg?branch=master
    :target: https://travis-ci.org/tanaylab/metacells
    :alt: Build Status

.. image:: https://readthedocs.org/projects/metacells/badge/?version=latest
    :target: https://metacells.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

*TODO*

The metacells package implements the improved metacell algorithm [2]_ for single-cell RNA sequencing
(scRNA-seq) data analysis within the `scipy https://www.scipy.org/` framework. The original metacell
algorithm [1]_ was implemented in R. The python package contains various algorithmic improvements
and is scalable for larger data sets (millions of cells).

Metacell Analysis
=================

Naively, scRNA_seq data is a set of cell profiles, where for each one, for each gene, we get a count
of the mRNA molecules that existed in the cell for that gene. This serves as an indicator of how
"expressed" or "active" the gene is.

As in any real world technology, the raw data may suffer from technical artifacts (counting the
molecules of two cells in one profile, counting the molecules from a ruptured cells, counting only
the molecules from the cell nucleus, etc.). This requires pruning the raw data to exclude such
artifacts.

The current technology scRNA-seq data is also very sparse (typically <<10% the RNA molecules are
counted). This introduces large sampling variance on top of the original signal, which itself
contains significant inherent biological noise.

Analyzing scRNA-seq data therefore requires processing the profiles in bulk. Classically, this has
been done by directly clustering the cells using various methods.

In contrast, the metacell approach groups together profiles of the "same" biological state into
groups of cells of the "same" biological state, with the *minimal* number of profiles needed for
computing robust statistics (in particular, mean gene expression). Each such group is a single
"metacell".

By summing profiles of cells of the "same" state together, each metacell greatly reduces the
sampling variance, and provides a more robust estimation of the transcription state. Note a metacell
is *not* a cell type (multiple metacells may belong to the same "type", or even have the "same"
state, if the data sufficiently over-samples this state). Also, a metacell is *not* a parametric
model of the cell state. It is merely a more robust description of some cell state.

The metacells should therefore be further analyzed as if they were cells, using additional methods
to classify cell types, detect cell trajectories and/or lineage, build parametric models for cell
behavior, etc. Using metacells as input for such analysis techniques should benefit both from the
more robust, less noisy input; and also from the (~100-fold) reduction in the number of cells to
analyze when dealing with large data (e.g. analyzing millions of individual cells).

Usage
-----

**TODO**

Installation
============

In short: ``pip install metacells``.

**TODO**: Alternatively, use `Anaconda <https://www.anaconda.com/distribution/>`_ and get the conda
packages from the `conda-forge channel <https://anaconda.org/conda-forge/metacells>`_, which
supports both Unix, Mac OS and Windows.

For Unix like systems it is possible to install from source. For Windows this is overly complicated,
and you are recommended to use the binary wheels. Also, this package is meant to be used in the
context of the `scipy <http://scipy.org>`_ framework.

Make sure you have all necessary tools for compilation. In Ubuntu this can be installed using ``sudo
apt-get install build-essential``, please refer to the documentation for your specific system.  Make
sure that not only ``gcc`` is installed, but also ``g++``, as the ``metacells`` package is
programmed in ``C++``.

You can check if all went well by running a variety of tests using ``pytest`` or ``tox``.

References
==========

Please cite the references appropriately in case they are used.

.. [1] Baran, Y., Bercovich, A., Sebe-Pedros, A. et al. MetaCell: analysis of single-cell RNA-seq
   data using K-nn graph partitions. Genome Biol 20, 206 (2019).
   `10.1186/s13059-019-1812-2 <https://doi.org/10.1186/s13059-019-1812-2>`_

.. [2] *TODO*.

License (MIT)
=============

Copyright © 2020 Weizmann Institute of Science

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

