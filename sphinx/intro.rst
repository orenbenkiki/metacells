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
