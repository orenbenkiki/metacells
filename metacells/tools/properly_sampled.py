'''
Properly Sampled
----------------
'''

import logging
from typing import Optional

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
from anndata import AnnData

import metacells.parameters as pr
import metacells.preprocessing as pp
import metacells.utilities as ut

__all__ = [
    'find_properly_sampled_cells',
    'find_properly_sampled_genes',
]


LOG = logging.getLogger(__name__)


@ut.timed_call()
@ut.expand_doc()
def find_properly_sampled_cells(
    adata: AnnData,
    of: Optional[str] = None,
    *,
    min_cell_total: Optional[int] = pr.properly_sampled_min_cell_total,
    max_cell_total: Optional[int] = pr.properly_sampled_max_cell_total,
    inplace: bool = True,
    intermediate: bool = True,
) -> Optional[ut.PandasSeries]:
    '''
    Detect cells with a "proper" amount ``of`` some data sampled (by default, the focus).

    Due to both technical effects and natural variance between cells, the total number of UMIs
    varies from cell to cell. We often would like to work on cells that contain a sufficient number
    of UMIs for meaningful analysis; we sometimes also wish to exclude cells which have "too many"
    UMIs.

    **Input**

    A :py:func:`metacells.utilities.annotation.setup` annotated ``adata``, where the observations
    are cells and the variables are genes.

    **Returns**

    Observation (Cell) Annotations
        ``properly_sampled_cell``
            A boolean mask indicating whether each cell has a "proper" amount of samples.

    If ``inplace`` (default: {inplace}), this is written to the data, and the function returns
    ``None``. Otherwise this is returned as a pandas series (indexed by the observation names).

    If ``intermediate`` (default: {intermediate}), keep all all the intermediate data (e.g. sums)
    for future reuse. Otherwise, discard it.

    **Computation Parameters**

    1. Exclude all cells whose total data is less than the ``min_cell_total`` (default:
       {min_cell_total}), unless it is ``None``

    2. Exclude all cells whose total data is more than the ``max_cell_total`` (default:
       {max_cell_total}), unless it is ``None``
    '''
    of, level = ut.log_operation(LOG, adata, 'find_properly_sampled_cells', of)

    with ut.focus_on(ut.get_vo_data, adata, of, intermediate=intermediate):
        total_of_cells = pp.get_per_obs(adata, ut.sum_per).proper
        cells_mask = np.full(adata.n_obs, True, dtype='bool')

        if min_cell_total is not None:
            LOG.debug('  min_cell_total: %s', min_cell_total)
            cells_mask = cells_mask & (total_of_cells >= min_cell_total)

        if max_cell_total is not None:
            LOG.debug('  max_cell_total: %s', max_cell_total)
            cells_mask = \
                cells_mask & (total_of_cells <= max_cell_total)

    if inplace:
        ut.set_o_data(adata, 'properly_sampled_cell', cells_mask)
        return None

    ut.log_mask(LOG, level, 'properly_sampled_cell', cells_mask)

    return pd.Series(cells_mask, index=adata.obs_names)


@ut.timed_call()
@ut.expand_doc()
def find_properly_sampled_genes(
    adata: AnnData,
    of: Optional[str] = None,
    *,
    min_gene_total: int = pr.properly_sampled_min_gene_total,
    inplace: bool = True,
    intermediate: bool = True,
) -> Optional[ut.PandasSeries]:
    '''
    Detect genes with a "proper" amount ``of`` some data samples (by default, the focus).

    Due to both technical effects and natural variance between genes, the expression of genes varies
    greatly between cells. This is exactly the information we are trying to analyze. We often would
    like to work on genes that have a sufficient level of expression for meaningful analysis.
    Specifically, it doesn't make sense to analyze genes that have zero expression in all the cells.

    .. todo::

        Provide additional optional criteria for "properly sampled genes"?

    **Input**

    A :py:func:`metacells.utilities.annotation.setup` annotated ``adata``, where the observations
    are cells and the variables are genes.

    **Returns**

    Variable (Gene) Annotations
        ``properly_sampled_gene``
            A boolean mask indicating whether each gene has a "proper" number of samples.

    If ``inplace`` (default: {inplace}), this is written to the data and the function returns
    ``None``. Otherwise this is returned as a pandas series (indexed by the variable names).

    If ``intermediate`` (default: {intermediate}), keep all all the intermediate data (e.g. sums)
    for future reuse. Otherwise, discard it.

    **Computation Parameters**

    1. Exclude all genes whose total data is less than the ``min_gene_total`` (default:
       {min_gene_total}).
    '''
    of, level = ut.log_operation(LOG, adata, 'find_properly_sampled_genes', of)

    with ut.focus_on(ut.get_vo_data, adata, of, intermediate=intermediate):
        total_of_genes = pp.get_per_var(adata, ut.sum_per).proper
        LOG.debug('  min_gene_total: %s', min_gene_total)
        genes_mask = total_of_genes >= min_gene_total

    if inplace:
        ut.set_v_data(adata, 'properly_sampled_gene', genes_mask)
        return None

    ut.log_mask(LOG, level, 'properly_sampled_gene', genes_mask)

    return pd.Series(genes_mask, index=adata.obs_names)
