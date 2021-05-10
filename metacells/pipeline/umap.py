'''
UMAP
----

It is useful to project the metacells data to a 2D scatter plot (where each point is a metacell).
The steps provided here are expected to yield a reasonable such projection, though as always you
might need to tweak the parameters or even the overall flow for specific data sets.
'''

from typing import Union

import numpy as np
import scipy.sparse as sparse  # type: ignore
from anndata import AnnData

import metacells.parameters as pr
import metacells.tools as tl
import metacells.utilities as ut

__all__ = [
    'compute_umap_by_features',
]


@ut.logged()
@ut.timed_call()
@ut.expand_doc()
def compute_umap_by_features(
    adata: AnnData,
    what: Union[str, ut.Matrix] = '__x__',
    *,
    max_top_feature_genes: int = pr.max_top_feature_genes,
    similarity_value_normalization: float = pr.umap_similarity_value_normalization,
    similarity_log_data: bool = pr.umap_similarity_log_data,
    similarity_method: str = pr.umap_similarity_method,
    logistics_location: float = pr.logistics_location,
    logistics_scale: float = pr.logistics_scale,
    k: int = pr.umap_k,
    min_dist: float = pr.umap_min_dist,
    spread: float = pr.umap_spread,
    random_seed: int = pr.random_seed,
) -> None:
    '''
    Compute a UMAP projection of the (meta)cells.

    **Input**

    Annotated ``adata`` where each observation is a metacells and the variables are genes,
    are genes, where ``what`` is a per-variable-per-observation matrix or the name of a
    per-variable-per-observation annotation containing such a matrix.

    **Returns**

    Sets the following annotations in ``adata``:

    Variable (Gene) Annotations
        ``top_feature_gene``
            A boolean mask of the top feature genes used to compute similarities between the
            metacells.

    Observation (Metacell) Annotations
        ``umap_x``, ``umap_y``
            The X and Y coordinates of each metacell in the UMAP projection.

    **Computation Parameters**

    1. Invoke :py:func:`metacells.tools.high.find_top_feature_genes` using ``max_top_feature_genes``
       (default: {max_top_feature_genes}) to pick the feature genes to use to compute similarities
       between the metacells.

    2. Compute the fractions of each gene in each cell, and add the
       ``similarity_value_normalization`` (default: {similarity_value_normalization}) to
       it.

    3. If ``similarity_log_data`` (default: {similarity_log_data}), invoke the
       :py:func:`metacells.utilities.computation.log_data` function to compute the log (base 2) of
       the data.

    4. Invoke :py:func:`metacells.tools.similarity.compute_obs_obs_similarity` using
       ``similarity_method`` (default: {similarity_method}), ``logistics_location`` (default:
       {logistics_location}) and ``logistics_scale`` (default: {logistics_scale}) and convert this
       to distances.

    5. Invoke :py:func:`metacells.tools.layout.umap_by_distances` using these distances, ``k``
       (default: {k}), ``min_dist`` (default: {min_dist}), ``spread`` (default: {spread}) and the
       ``random_seed`` (default: {random_seed}). Note that if the seed is not zero, the result will
       be reproducible, but the computation will use only a single thread which will take longer to
       complete.
    '''
    tl.find_top_feature_genes(adata, max_genes=max_top_feature_genes)

    all_data = ut.get_vo_proper(adata, what, layout='row_major')
    all_fractions = ut.fraction_by(all_data, by='row')

    top_feature_genes_mask = ut.get_v_numpy(adata, 'top_feature_gene')

    top_feature_genes_fractions = all_fractions[:, top_feature_genes_mask]
    top_feature_genes_fractions = \
        ut.to_layout(top_feature_genes_fractions, layout='row_major')

    top_feature_genes_fractions += similarity_value_normalization

    if similarity_log_data:
        top_feature_genes_fractions = \
            ut.log_data(top_feature_genes_fractions, base=2)

    tdata = ut.slice(adata, vars=top_feature_genes_mask)
    similarities = tl.compute_obs_obs_similarity(tdata,
                                                 top_feature_genes_fractions,
                                                 method=similarity_method,
                                                 logistics_location=logistics_location,
                                                 logistics_scale=logistics_scale,
                                                 inplace=False)
    assert similarities is not None

    distances = ut.to_numpy_matrix(similarities)
    distances *= -1
    distances += 1
    np.fill_diagonal(distances, 0.0)
    distances = sparse.csr_matrix(distances)

    tl.umap_by_distances(adata, distances, k=k, min_dist=min_dist,
                         spread=spread, random_seed=random_seed)