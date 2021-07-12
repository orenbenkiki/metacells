'''
Feature
-------
'''

from re import Pattern
from typing import Collection, Optional, Union

from anndata import AnnData

import metacells.parameters as pr
import metacells.tools as tl
import metacells.utilities as ut

__all__ = [
    'extract_feature_data',
]


@ut.logged()
@ut.timed_call()
@ut.expand_doc()
def extract_feature_data(
    adata: AnnData,
    what: Union[str, ut.Matrix] = '__x__',
    *,
    name: Optional[str] = '.feature',
    downsample_min_samples: float = pr.feature_downsample_min_samples,
    downsample_min_cell_quantile: float = pr.feature_downsample_min_cell_quantile,
    downsample_max_cell_quantile: float = pr.feature_downsample_max_cell_quantile,
    min_gene_relative_variance: Optional[float] = pr.feature_min_gene_relative_variance,
    min_gene_total: Optional[int] = pr.feature_min_gene_total,
    min_gene_top3: Optional[int] = pr.feature_min_gene_top3,
    forbidden_gene_names: Optional[Collection[str]] = None,
    forbidden_gene_patterns: Optional[Collection[Union[str, Pattern]]] = None,
    random_seed: int = 0,
    top_level: bool = True,
) -> Optional[AnnData]:
    '''
    Extract a "feature" subset of ``what`` (default: {what} data, to compute metacells by.

    When computing metacells (or clustering cells in general), it makes sense to use a subset of the
    genes for computing cell-cell similarity, for both technical (e.g., too low an expression level)
    and biological (e.g., ignoring bookkeeping and cell cycle genes) reasons. The steps provided
    here are expected to be generically useful, but as always specific data sets may require custom
    feature selection steps on a case-by-case basis.

    **Input**

    A presumably "clean" Annotated ``adata``, where the observations are cells and the variables are
    genes, where ``what`` is a per-variable-per-observation matrix or the name of a
    per-variable-per-observation annotation containing such a matrix.

    **Returns**

    Returns annotated sliced data containing the "feature" subset of the original data. By default,
    the ``name`` of this data is {name}. If no features were selected, return ``None``.

    Also sets the following annotations in the full ``adata``:

    Variable (Gene) Annotations
        ``high_total_gene``
            A boolean mask of genes with "high" expression level.

        ``high_relative_variance_gene``
            A boolean mask of genes with "high" normalized variance, relative to other genes with a
            similar expression level.

        ``forbidden_gene``
            A boolean mask of genes which are forbidden from being chosen as "feature" genes based
            on their name.

        ``feature_gene``
            A boolean mask of the "feature" genes.

    **Computation Parameters**

    1. Invoke :py:func:`metacells.tools.downsample.downsample_cells` to downsample the cells to the
       same total number of UMIs, using the ``downsample_min_samples`` (default:
       {downsample_min_samples}), ``downsample_min_cell_quantile`` (default:
       {downsample_min_cell_quantile}), ``downsample_max_cell_quantile`` (default:
       {downsample_max_cell_quantile}) and the ``random_seed`` (default: {random_seed}).

    2. Invoke :py:func:`metacells.tools.high.find_high_total_genes` to select high-expression
       feature genes (based on the downsampled data), using ``min_gene_total``.

    3. Invoke :py:func:`metacells.tools.high.find_high_relative_variance_genes` to select
       high-variance feature genes (based on the downsampled data), using
       ``min_gene_relative_variance``.

    4. Invoke :py:func:`metacells.tools.named.find_named_genes` to forbid genes from being used as
       feature genes, based on their name. using the ``forbidden_gene_names`` (default:
       {forbidden_gene_names}) and ``forbidden_gene_patterns`` (default: {forbidden_gene_patterns}).
       This is stored in an intermediate per-variable (gene) ``forbidden_genes`` boolean mask.

    5. Invoke :py:func:`metacells.tools.filter.filter_data` to slice just the selected
       "feature" genes using the ``name`` (default: {name}).
    '''
    tl.downsample_cells(adata, what,
                        downsample_min_samples=downsample_min_samples,
                        downsample_min_cell_quantile=downsample_min_cell_quantile,
                        downsample_max_cell_quantile=downsample_max_cell_quantile,
                        random_seed=random_seed)

    var_masks = []

    if min_gene_top3 is not None:
        var_masks.append('high_top3_gene')
        tl.find_high_topN_genes(adata, 'downsampled', topN=3,
                                min_gene_topN=min_gene_top3)

    if min_gene_total is not None:
        var_masks.append('high_total_gene')
        tl.find_high_total_genes(adata, 'downsampled',
                                 min_gene_total=min_gene_total)

    if min_gene_relative_variance is not None:
        var_masks.append('high_relative_variance_gene')
        tl.find_high_relative_variance_genes(adata, 'downsampled',
                                             min_gene_relative_variance=min_gene_relative_variance)

    if forbidden_gene_names is not None \
            or forbidden_gene_patterns is not None:
        var_masks.append('~forbidden_gene')
        tl.find_named_genes(adata,
                            to='forbidden_gene',
                            names=forbidden_gene_names,
                            patterns=forbidden_gene_patterns)

    results = tl.filter_data(adata, name=name, top_level=top_level,
                             track_var='full_gene_index',
                             mask_var='feature_gene',
                             var_masks=var_masks)

    if results is None:
        raise ValueError('Empty feature data, giving up')

    return results[0]
