'''
Mask genes by their name.
'''

from re import Pattern
from typing import Optional, Union

import pandas as pd  # type: ignore
from anndata import AnnData

import metacells.utilities as ut

__all__ = [
    'find_named_genes',
]


@ut.timed_call()
@ut.expand_doc()
def find_named_genes(
    adata: AnnData,
    pattern: Union[str, Pattern],
    *,
    name: Optional[str] = None,
    invert: bool = False,
) -> Optional[ut.PandasSeries]:
    '''
    Find genes by their name.

    This creates a mask of all the genes whose name matches the ``pattern``. If ``invert`` (default:
    {invert}), invert the resulting mask.

    If ``name`` (default: {name}) is specified, this is stored as a per-variable (gene) annotation
    with that name, and returns ``None``. This is useful to fill gene masks such as
    ``excluded_genes`` (genes which should be excluded from the rest of the processing) and
    ``forbidden_genes`` (genes which must not be chosen as feature genes).

    Otherwise, it returns it as a pandas series (indexed by the variable, that is gene, names).
    '''
    genes_mask = ut.regex_matches_mask(pattern, adata.var_names)
    if invert:
        genes_mask = ~genes_mask

    if name is None:
        return pd.Series(genes_mask, index=adata.var_names)

    adata.var[name] = genes_mask
    ut.safe_slicing_data(name, ut.ALWAYS_SAFE)
    return None
