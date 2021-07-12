Metacell Analysis
=================

**DANGER WILL ROBINSON:** This is a pre-alpha version under active development, so please do not use
this (yet). We are working on a 1st release, hopefully in the next few weeks, in combination with
submitting a paper.

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

Installation
============

In short: ``pip install metacells``. If you do not have ``sudo`` privileges, you might need to ``pip
install --user metacells``. Note that ``metacells`` requires many "heavy" dependencies, most notably
``numpy``, ``pandas``, ``scipy``, ``scanpy``, which ``pip`` should automatically install for you.

The metacells package contains extensions written in C++, which means installing it requires
compilation. In Ubuntu the necessary tools can be installed using ``sudo apt-get install
build-essential`` - please refer to the documentation for your specific system. Make sure that not
only ``gcc`` is installed, but also ``g++`` to allow for C++ compilation.

You can check if all went well by running a variety of tests using ``pytest`` or ``tox``.

Vignettes
=========

The `generated documentation <https://metacells.readthedocs.io/en/latest>`_
contains the following vignettes:
`Basic Metacells Vignette <https://metacells.readthedocs.io/en/latest/Metacells_Vignette.html>`_,
`Manual Analysis Vignette <https://metacells.readthedocs.io/en/latest/Manual_Analysis.html>`_,
and
`Seurat Analysis Vignette <https://metacells.readthedocs.io/en/latest/Seurat_Analysis.html>`_.

You can also access their very latest version in the `Github repository
<https://github.com/tanaylab/metacells>`_ in:
`Basic Metacells Vignette <https://github.com/tanaylab/metacells/blob/master/sphinx/Metacells_Vignette.rst>`_,
`Manual Analysis Vignette <https://github.com/tanaylab/metacells/blob/master/sphinx/Manual_Analysis.rst>`_,
and
`Seurat Analysis Vignette <https://github.com/tanaylab/metacells/blob/master/sphinx/Seurat_Analysis.rst>`_.
