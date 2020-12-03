'''
Distincts
---------
'''

import logging
from typing import Optional, Tuple, Union

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
from anndata import AnnData

import metacells.extensions as xt  # type: ignore
import metacells.parameters as pr
import metacells.preprocessing as pp
import metacells.utilities as ut

__all__ = [
    'compute_distinct_folds',
    'find_distinct_genes',
    'compute_subset_distinct_genes',
]


LOG = logging.getLogger(__name__)


@ut.timed_call()
@ut.expand_doc()
def compute_distinct_folds(
    adata: AnnData,
    *,
    of: Optional[str] = None,
    normalization: float = 0,
    inplace: bool = True,
    intermediate: bool = True,
) -> Optional[ut.PandasFrame]:
    '''
    Compute for each observation (cell) and each variable (gene) how much is the value
    different from the overall population.

    **Input**

    A :py:func:`metacells.utilities.annotation.setup` annotated ``adata``, where the observations
    are (mata)cells and the variables are genes.

    **Returns**

    Per-Observation-Per-Variable (Cell-Gene) Annotations:
        ``distinct_ratio``
            For each gene in each cell, the log based 2 of the ratio between the fraction of the
            gene in the cell and the fraction of the gene in the overall population (sum of cells).

    If ``inplace`` (default: {inplace}), this is written to the data, and the function returns
    ``None``. Otherwise this is returned as a pandas frame (indexed by the observation and
    distinct gene rank).

    If ``intermediate`` (default: {intermediate}), keep all all the intermediate data (e.g. sums)
    for future reuse. Otherwise, discard it.

    **Computation Parameters**

    1. Compute, for each gene, the fraction of the gene's values out of the total sum of the values
       (that is, the mean fraction of the gene's expression in the population).

    2. Compute, for each cell, for each gene, the fraction of the gene's value out of the sum
       of the values in the cell (that is, the fraction of the gene's expression in the cell).

    3. Divide the two to the distinct ratio (that is, how much the gene's expression in the cell is
       different from the overall population), first adding the ``normalization`` (default:
       {normalization}) to both.

    4. Compute the log (base 2) of the result and use it as the fold factor.
    '''
    of, _ = ut.log_operation(LOG, adata, 'compute_distinct_folds', of)

    with ut.focus_on(ut.get_vo_data, adata, of, layout='row_major',
                     intermediate=intermediate) as data:
        total_umis_of_cells = pp.get_per_obs(adata, ut.sum_per).proper
        total_umis_of_genes = pp.get_per_var(adata, ut.sum_per).proper

        fractions_of_genes_in_data = total_umis_of_genes / \
            np.sum(total_umis_of_genes)
        fractions_of_genes_in_data += normalization

        total_umis_of_cells = np.copy(total_umis_of_cells)
        total_umis_of_cells[total_umis_of_cells == 0] = 1
        fraction_of_genes_in_cells = data / total_umis_of_cells[:, np.newaxis]
        fraction_of_genes_in_cells += normalization

        zeros_mask = fractions_of_genes_in_data <= 0
        fractions_of_genes_in_data[zeros_mask] = -1
        fraction_of_genes_in_cells[:, zeros_mask] = -1

        ratio_of_genes_in_cells = fraction_of_genes_in_cells
        ratio_of_genes_in_cells /= fractions_of_genes_in_data
        assert np.min(np.min(ratio_of_genes_in_cells)) > 0

        fold_of_genes_in_cells = \
            np.log2(ratio_of_genes_in_cells, out=ratio_of_genes_in_cells)

    if inplace:
        ut.set_vo_data(adata, 'distinct_fold', fold_of_genes_in_cells)
        return None

    return pd.DataFrame(fold_of_genes_in_cells, index=adata.obs_names, columns=adata.var_names)


@ut.timed_call()
@ut.expand_doc()
def find_distinct_genes(
    adata: AnnData,
    *,
    distinct_fold: str = 'distinct_fold',
    distinct_genes_count: int = pr.distinct_genes_count,
    inplace: bool = True,
) -> Optional[Tuple[ut.PandasFrame, ut.PandasFrame]]:
    '''
    Find for each observation (cell) the genes in which it is most distinct from the general
    population. This is typically applied to the metacells data rather than to the cells data.

    **Input**

    A :py:func:`metacells.utilities.annotation.setup` annotated ``adata``, where the observations
    are (mata)cells and the variables are genes, including a per-observation-per-variable annotated
    data named ``distinct_fold``, e.g. as computed by :py:func:`compute_distinct_folds`.

    **Returns**

    Observation-Any (Cell) Annotations
        ``cell_distinct_gene_indices``
            For each cell, the indices of its top ``distinct_genes_count`` genes.
        ``cell_distinct_gene_folds``
            For each cell, the fold factor of its top ``distinct_genes_count``.

    If ``inplace`` (default: {inplace}), this is written to the data, and the function returns
    ``None``. Otherwise this is returned as two pandas frames (indexed by the observation and
    distinct gene rank).

    **Computation Parameters**

    1. Fetch the previously computed per-observation-per-variable `distinct_fold`` annotation
       (default: {distinct_fold}).

    2. Keep the ``distinct_genes_count`` (default: {distinct_genes_count}) top fold factors.
    '''
    ut.log_operation(LOG, adata, 'find_distinct_genes')
    assert 0 < distinct_genes_count < adata.n_vars

    distinct_gene_indices = \
        np.empty((adata.n_obs, distinct_genes_count), dtype='int32')
    distinct_gene_folds = \
        np.empty((adata.n_obs, distinct_genes_count), dtype='float32')

    fold_in_cells = ut.get_vo_data(adata, distinct_fold, layout='row_major')
    xt.top_distinct(distinct_gene_indices, distinct_gene_folds,
                    fold_in_cells, False)

    if inplace:
        ut.set_oa_data(adata, 'cell_distinct_gene_indices',
                       distinct_gene_indices)
        ut.set_oa_data(adata, 'cell_distinct_gene_folds', distinct_gene_folds)
        return None

    return pd.DataFrame(distinct_gene_indices, index=adata.obs_names), \
        pd.DataFrame(distinct_gene_folds, index=adata.obs_names)


@ut.timed_call()
@ut.expand_doc()
def compute_subset_distinct_genes(
    adata: AnnData,
    *,
    to: Optional[str] = None,
    distinct_fold: str = 'distinct_fold',
    subset: Union[str, ut.DenseVector],
    included_fold_quantile: float = 0,
    excluded_fold_quantile: float = 1,
    inplace: bool = True,
) -> Optional[ut.PandasSeries]:
    '''
    Given a subset of the observations (cells), compute for each gene how distinct it is
    in the subset compared to the overall population.

    This is the difference between the gene's fold factor within the subset to its fold factor in
    the complement set. For example, for some gene, if the subset contains observations which are
    distinct from the population, but the complement set also contains observations which are
    distinct from the population, then this gene doesn't really distinguish the subset from its
    complement set. The quantile parameters control how strong this effect is; by default, we take
    the minimal fold factor and subtract the maximal fold factor in the complement set, which might
    be too strict.

    **Input**

    A :py:func:`metacells.utilities.annotation.setup` annotated ``adata``, where the observations
    are (mata)cells and the variables are genes, including a per-observation-per-variable annotated
    data named ``distinct_fold``, e.g. as computed by :py:func:`compute_distinct_folds`.

    **Returns**

    Variable (Gene) Annotations
        ``to``
            Store the fold factor of the gene in the subset as opposed to the overall population.

    If ``inplace`` (default: {inplace}), this is written to the data, and the function returns
    ``None``. This requires ``to`` to be specified. Otherwise this is returned as a pandas series
    (indexed by the gene names).

    **Computation Parameters**

    1. Fetch the previously computed per-observation-per-variable `distinct_fold`` annotation
       (default: {distinct_fold}).

    2. Split the data to the ``subset`` (either the name of a boolean mask, per-observation
       annotation, or a boolean mask vector, or a vector of observation indices) and its
       complement set.

    3. For each gene, compute the ``included_fold_quantile`` (default: {included_fold_quantile})
       of the fold factors of the gene in the cells in the subset and the ``excluded_fold_quantile``
       (default: {excluded_fold_quantile}) of the fold factors of the genes in the cells in the
       complement set.

    4. The difference between the two is the excess fold factor for the gene for the subset.
    '''
    ut.log_operation(LOG, adata, 'find_subset_distinct_genes')

    if isinstance(subset, str):
        subset = ut.to_proper_vector(ut.get_o_data(adata, subset))

    if subset.dtype != 'bool':
        mask: ut.DenseVector = np.full(False, adata.n_obs)
        mask[subset] = True
        subset = mask

    complement = ~subset

    fold_in_cells = ut.get_vo_data(adata, distinct_fold, layout='row_major')
    fold_in_subset = fold_in_cells[subset, :]
    fold_in_complement = fold_in_cells[complement, :]
    quantile_in_subset = \
        np.quantile(fold_in_subset, included_fold_quantile, axis=0)
    quantile_in_complement = \
        np.quantile(fold_in_complement, excluded_fold_quantile, axis=0)
    subset_distinct_fold = quantile_in_subset - quantile_in_complement

    if inplace:
        assert to is not None
        ut.set_v_data(adata, to, subset_distinct_fold)
        return None

    return pd.Series(subset_distinct_fold, index=adata.var_names)
