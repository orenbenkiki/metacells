'''
Utilities for dealing with data annotations.

This tries to deal safely with data layers, observation, variable and unstructured annotations,
especially in the presence to slicing the data.

.. todo::

    The :py:mod:`metacells.utilities.annotation` module should also deal with multi-dimensional
    ``obsm`` and ``varm`` annotations.

'''

from typing import Any, Callable, Collection, Dict, Optional, Set, Tuple, Union
from warnings import warn

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import scipy.sparse as sparse  # type: ignore
from anndata import AnnData  # type: ignore
from readerwriterlock import rwlock

import metacells.utilities.computation as utc
import metacells.utilities.documentation as utd
import metacells.utilities.timing as timed

__all__ = [
    'slice',
    'slicing_data_layer',
    'slicing_obs_annotation',
    'slicing_var_annotation',
    'slicing_uns_annotation',

    'annotate_as_base',

    'set_x_layer',

    'get_data_layer',
    'del_data_layer',
    'get_log_data',

    'get_annotation_of_obs',
    'get_sum_per_obs',
    'get_nnz_per_obs',
    'get_max_per_obs',

    'get_annotation_of_var',
    'get_sum_per_var',
    'get_nnz_per_var',
    'get_max_per_var',

    'get_annotation_of_data',
]


LOCK = rwlock.RWLockRead()

WARN_WHEN_SLICING_UNKNOWN_ANNOTATIONS: bool = True

SAFE_SLICING_OBS_ANNOTATIONS: Dict[str, Tuple[bool, bool]] = {}
SAFE_SLICING_VAR_ANNOTATIONS: Dict[str, Tuple[bool, bool]] = {}
SAFE_SLICING_UNS_ANNOTATIONS: Dict[str, Tuple[bool, bool]] = {}
SAFE_SLICING_DATA_LAYERS: Dict[str, Tuple[bool, bool]] = {}


def slicing_data_layer(
    name: str,
    *,
    preserve_when_slicing_obs: bool = False,
    preserve_when_slicing_var: bool = False,
) -> None:
    '''
    Specify whether the named data layer should be preserved when slicing either the observations
    (cells) or variables (genes).

    When the data is sliced, a data layer may become invalidated and should not be preserved. For
    example, the ``sum_UMIs`` of genes is invalidated when specifying a subset of the cells.
    '''
    with LOCK.gen_wlock():
        SAFE_SLICING_DATA_LAYERS[name] = \
            (preserve_when_slicing_obs, preserve_when_slicing_var)


def _slicing_known_data_layers() -> None:
    slicing_data_layer('UMIs',
                       preserve_when_slicing_obs=True,
                       preserve_when_slicing_var=True)
    slicing_data_layer('log_UMIs',
                       preserve_when_slicing_obs=True,
                       preserve_when_slicing_var=True)
    slicing_data_layer('downsampled_UMIs',
                       preserve_when_slicing_obs=True,
                       preserve_when_slicing_var=False)


def slicing_obs_annotation(
    name: str,
    *,
    preserve_when_slicing_obs: bool = False,
    preserve_when_slicing_var: bool = False,
) -> None:
    '''
    Specify whether the named per-observation (cell) annotation should be preserved when slicing
    either the observations (cells) or variables (genes).
    '''
    SAFE_SLICING_OBS_ANNOTATIONS[name] = \
        (preserve_when_slicing_obs, preserve_when_slicing_var)


def _slicing_known_obs_annotations() -> None:
    slicing_obs_annotation('base_index',
                           preserve_when_slicing_obs=True,
                           preserve_when_slicing_var=True)
    slicing_obs_annotation('sum_UMIs',
                           preserve_when_slicing_obs=True,
                           preserve_when_slicing_var=False)
    slicing_obs_annotation('nnz_UMIs',
                           preserve_when_slicing_obs=True,
                           preserve_when_slicing_var=False)
    slicing_obs_annotation('max_UMIs',
                           preserve_when_slicing_obs=True,
                           preserve_when_slicing_var=False)


def slicing_var_annotation(
    name: str,
    *,
    preserve_when_slicing_obs: bool = False,
    preserve_when_slicing_var: bool = False,
) -> None:
    '''
    Specify whether the named per-variable (gene) annotation should be preserved when slicing
    either the observations (cells) or variables (genes).
    '''
    SAFE_SLICING_VAR_ANNOTATIONS[name] = \
        (preserve_when_slicing_obs, preserve_when_slicing_var)


def _slicing_known_var_annotations() -> None:
    slicing_var_annotation('gene_ids',
                           preserve_when_slicing_obs=True,
                           preserve_when_slicing_var=True)
    slicing_var_annotation('base_index',
                           preserve_when_slicing_obs=True,
                           preserve_when_slicing_var=True)
    slicing_var_annotation('sum_UMIs',
                           preserve_when_slicing_obs=False,
                           preserve_when_slicing_var=True)
    slicing_obs_annotation('nnz_UMIs',
                           preserve_when_slicing_obs=False,
                           preserve_when_slicing_var=True)
    slicing_obs_annotation('max_UMIs',
                           preserve_when_slicing_obs=False,
                           preserve_when_slicing_var=True)


def slicing_uns_annotation(
    name: str,
    *,
    preserve_when_slicing_obs: bool = False,
    preserve_when_slicing_var: bool = False,
) -> None:
    '''
    Specify whether the named unstructured annotation should be preserved when slicing either the
    observations (cells) or variables (genes).
    '''
    with LOCK.gen_wlock():
        SAFE_SLICING_UNS_ANNOTATIONS[name] = \
            (preserve_when_slicing_obs, preserve_when_slicing_var)


def _slicing_known_uns_annotations() -> None:
    slicing_uns_annotation('x_layer',
                           preserve_when_slicing_obs=True,
                           preserve_when_slicing_var=True)


def _slicing_known_annotations() -> None:
    _slicing_known_data_layers()
    _slicing_known_obs_annotations()
    _slicing_known_var_annotations()
    _slicing_known_uns_annotations()


_slicing_known_annotations()


@timed.call()
def slice(  # pylint: disable=redefined-builtin
    adata: AnnData,
    *,
    cells: Optional[Collection] = None,
    genes: Optional[Collection] = None,
    invalidated_annotations_prefix: Optional[str] = None,
) -> AnnData:
    '''
    Return new annotated data which includes a subset of the full ``adata``.

    If ``cells`` and/or ``genes`` are specified, they should include either a boolean
    mask or a collection of indices to include in the data slice. In the case of an indices array,
    it is assumed the indices are unique and sorted, that is that their effect is similar to a mask.

    In general, annotations data might become invalid when slicing (e.g., the ``sum_UMIs`` of genes
    is invalidated when specifying a subset of the cells). Therefore, annotations will be removed
    from the result unless they were explicitly marked as preserved using
    :py:func:`slicing_obs_annotation`, :py:func:`slicing_var_annotation`, and/or
    :py:func:`slicing_uns_annotation`.

    If ``invalidated_annotations_prefix`` is specified, then unpreserved annotations will not be
    removed; instead they will be renamed with the addition of the prefix.

    .. note::

        Setting the prefix to the empty string would preserves all the annotations including the
        invalidated ones. As this is unsafe, it will trigger a run-time exception. If you wish to
        perform such an unsafe slicing operation, invoke the built-in ``adata[..., ...]``.
    '''
    assert invalidated_annotations_prefix != ''

    if cells is None:
        cells = range(adata.n_obs)
    if genes is None:
        genes = range(adata.n_vars)

    will_slice_obs = len(cells) != adata.n_obs
    will_slice_var = len(genes) != adata.n_vars

    saved_data_layers = \
        _save_data_layers(adata, will_slice_obs, will_slice_var,
                          invalidated_annotations_prefix is not None)

    with timed.step('.data_layers'):
        bdata = adata[cells, genes]

    did_slice_obs = bdata.n_obs != adata.n_obs
    did_slice_var = bdata.n_vars != adata.n_vars

    assert did_slice_obs == will_slice_obs
    assert did_slice_var == will_slice_var

    if did_slice_obs or did_slice_var:
        _filter_layers(bdata, did_slice_obs, did_slice_var,
                       invalidated_annotations_prefix)
        _filter_annotations('obs', bdata.obs, did_slice_obs, did_slice_var,
                            SAFE_SLICING_OBS_ANNOTATIONS, invalidated_annotations_prefix)
        _filter_annotations('var', bdata.var, did_slice_obs, did_slice_var,
                            SAFE_SLICING_VAR_ANNOTATIONS, invalidated_annotations_prefix)
        _filter_annotations('uns', bdata.uns, did_slice_obs, did_slice_var,
                            SAFE_SLICING_UNS_ANNOTATIONS, invalidated_annotations_prefix)

    adata.layers.update(saved_data_layers)

    return bdata


def _save_data_layers(
    adata: AnnData,
    will_slice_obs: bool,
    will_slice_var: bool,
    will_prefix_invalidated: bool
) -> Dict[str, Any]:
    saved_data_layers: Dict[str, Any] = {}
    delete_data_layers: Set[str] = set()

    will_slice_only_obs = will_slice_obs and not will_slice_var
    will_slice_only_var = will_slice_var and not will_slice_obs

    for name, data in adata.layers.items():
        if name == adata.uns['x_layer']:
            delete_data_layers.add(name)

        else:
            action = \
                _analyze_data_layer(name, will_slice_obs,
                                    will_slice_var, will_prefix_invalidated)
            if action == 'discard':
                delete_data_layers.add(name)
                saved_data_layers[name] = data
                continue

            if name.startswith('__'):
                continue

        for will_slice_only, per in [(will_slice_only_obs, '__per_obs__'),
                                     (will_slice_only_var, '__per_var__')]:
            if not will_slice_only:
                continue

            per_name = per + name
            per_data = adata.layers.get(per_name)
            if per_data is None:
                break

            saved_data_layers[name] = data
            adata.layers[name] = per_data
            delete_data_layers.add(per_name)

            if name == adata.uns['x_layer']:
                with timed.step('.swap_x'):
                    adata.X = per_data

    for name in delete_data_layers:
        del adata.layers[name]

    return saved_data_layers


def _filter_layers(
    bdata: AnnData,
    did_slice_obs: bool,
    did_slice_var: bool,
    invalidated_annotations_prefix: Optional[str],
) -> None:
    x_layer = bdata.uns['x_layer']
    assert x_layer not in bdata.layers.keys()
    with timed.step('.x_layer'):
        bdata.layers[x_layer] = bdata.X

    for name in bdata.layers.keys():
        action = \
            _analyze_data_layer(name, did_slice_obs, did_slice_var,
                                invalidated_annotations_prefix is not None)

        assert action != 'discard'

        if action == 'preserve':
            data = bdata.layers[name]
            utc.freeze(data)
            continue

        data = bdata.layers[name]
        utc.freeze(data)

        assert action == 'prefix'
        assert invalidated_annotations_prefix is not None

        new_name = invalidated_annotations_prefix + name
        assert new_name != name
        bdata.layers[new_name] = data
        _clone_slicing_annotation(SAFE_SLICING_DATA_LAYERS, name, new_name)


def _analyze_data_layer(
    name: str,
    do_slice_obs: bool,
    do_slice_var: bool,
    will_prefix_invalidated: bool,
) -> str:
    if name[:2] != '__':
        is_per = False
        is_per_obs = False
        is_per_var = False
        base_name = name
    else:
        is_per = True
        is_per_obs = name.startswith('__per_obs__')
        is_per_var = name.startswith('__per_var__')
        assert is_per_obs or is_per_var
        base_name = name[11:]

    preserve_when_slicing_obs, preserve_when_slicing_var = \
        _get_slicing_annotation('data', 'layer',
                                SAFE_SLICING_DATA_LAYERS, base_name)

    if is_per_var:
        preserve_when_slicing_obs = False
    if is_per_obs:
        preserve_when_slicing_var = False

    preserve = (not do_slice_obs or preserve_when_slicing_obs) \
        and (not do_slice_var or preserve_when_slicing_var)

    if preserve:
        return 'preserve'

    if is_per or not will_prefix_invalidated:
        return 'discard'

    return 'prefix'


def _filter_annotations(  # pylint: disable=too-many-locals
    kind: str,
    new_annotations: Union[pd.DataFrame, Dict[str, Any]],
    did_slice_obs: bool,
    did_slice_var: bool,
    slicing_annotations: Dict[str, Tuple[bool, bool]],
    invalidated_annotations_prefix: Optional[str],
) -> None:
    annotation_names = list(new_annotations.keys())
    for name in annotation_names:
        preserve_when_slicing_obs, preserve_when_slicing_var = \
            _get_slicing_annotation(kind, 'annotation',
                                    slicing_annotations, name)

        preserve = (not did_slice_obs or preserve_when_slicing_obs) \
            and (not did_slice_var or preserve_when_slicing_var)

        if preserve:
            data = new_annotations[name]
            if kind != 'uns' \
                    or isinstance(data, (np.ndarray, sparse.spmatrix,
                                         pd.Series, pd.DataFrame)):
                utc.freeze(data)
            continue

        data = new_annotations[name]
        if kind != 'uns' \
                or isinstance(data, (np.ndarray, sparse.spmatrix,
                                     pd.Series, pd.DataFrame)):
            utc.freeze(data)

        if invalidated_annotations_prefix is not None:
            new_name = invalidated_annotations_prefix + name
            assert new_name != name
            new_annotations[new_name] = data
            _clone_slicing_annotation(slicing_annotations, name, new_name)

        del new_annotations[name]


def _get_slicing_annotation(
    kind: str,
    what: str,
    slicing_annotations: Dict[str, Tuple[bool, bool]],
    name: str
) -> Tuple[bool, bool]:
    with LOCK.gen_rlock():
        slicing_annotation = slicing_annotations.get(name)

    if slicing_annotation is None:
        with LOCK.gen_wlock():
            slicing_annotation = slicing_annotations.get(name)
            if slicing_annotation is None:
                unknown_sliced_annotation = \
                    'Slicing an unknown {kind} {what}: {name}; ' \
                    'assuming it should not be preserved' \
                    .format(kind=kind, what=what, name=name)
                warn(unknown_sliced_annotation)
                slicing_annotation = (False, False)
                slicing_annotations[name] = slicing_annotation

    return slicing_annotation


def _clone_slicing_annotation(
    slicing_annotations: Dict[str, Tuple[bool, bool]],
    old_name: str,
    new_name: str,
) -> None:
    with LOCK.gen_rlock():
        if new_name in slicing_annotations:
            return

    with LOCK.gen_wlock():
        if new_name in slicing_annotations:
            return
        slicing_annotations[new_name] = slicing_annotations[old_name]


def annotate_as_base(adata: AnnData, *, name: str = 'base_index') -> None:
    '''
    Create ``name`` (by default, ``base_index``) per-observation (cells) and per-variable (genes)
    annotations, which will be preserved when creating a :py:func:`slice` of the data to easily
    refer back to the original full data.
    '''
    adata.obs[name] = np.arange(adata.n_obs)
    adata.var[name] = np.arange(adata.n_vars)


def set_x_layer(adata: AnnData, name: str) -> None:
    '''
    Set the layer name of the data that is contained in ``X``. This is stored in the unstructured
    annotation named ``x_layer``.

    This should be called after populating the ``X`` data for the first time (e.g., by importing the
    data). It lets rest of the code to know what kind of data it holds. All the other layer
    utilities ``assert`` this was done.

    .. note::

        When using the layer utilities, do not directly read or write the value of ``X``. Instead
        use :py:func:`get_data_layer`.

    .. todo::

        Better integration with ``AnnData`` would allow accessing and even setting ``X`` in addition
        to using :py:func:`get_data_layer`. Currently this isn't implemented since ``AnnData``
        provides its own magic for handling ``X`` which is incompatible with the layers approach.
    '''
    X = adata.X
    assert X is not None
    assert 'x_layer' not in adata.uns_keys()

    adata.uns['x_layer'] = name
    adata.layers[name] = X


@timed.call()
def get_data_layer(
    adata: AnnData,
    name: str,
    compute: Optional[Callable[[], utc.Vector]] = None,
    *,
    inplace: bool = False,
    per: Optional[str] = None,
) -> utc.Matrix:
    '''
    Lookup a per-obsevration-per-variable matrix in ``adata`` by its ``name``.

    If the layer does not exist, ``compute`` it. If no ``compute`` function was given, ``raise``.

    If ``inplace``, store the annotation in ``adata``.

    If ``per`` is specified, it must be one of ``obs`` (cells) or ``var`` (genes). This returns the
    data in a layout optimized for per-observation (row-major / csr) or per-variable (column-major
    / csc). If also ``inplace``, this is cached in an additional layer whose name is prefixed (e.g.
    ``__per_obs__UMIs`` for the ``UMIs`` layer).

    .. note::

        In general layer names that start with ``__`` are reserved and should not be explicitly
        used.

    .. note::

        The original data returned by ``compute`` is always preserved under its non-prefixed layer
        name.

    .. todo::

        A better implementation of :py:func:`get_data_layer` would be to cache the layout-specific
        data in a private variable, but we do not own the ``AnnData`` object.
    '''

    assert name[:2] != '__'

    data = adata.layers.get(name)
    if data is not None:
        if name != adata.uns['x_layer']:
            utc.freeze(data)
    else:
        if compute is None:
            raise RuntimeError('unavailable layer data: ' + name)
        with timed.step('.compute'):
            data = compute()
        assert data.shape == adata.shape

        if inplace:
            utc.freeze(data)
            adata.layers[name] = data

    if per is not None:
        assert per in ['var', 'obs']
        per_name = '__per_%s__%s' % (per, name)
        axis = ['var', 'obs'].index(per)

        per_data = adata.layers.get(per_name)
        if per_data is None:
            per_data = utc.to_layout(data, axis=axis)

            if inplace:
                utc.freeze(per_data)
                adata.layers[per_name] = per_data

        data = per_data

    return data


def del_data_layer(
    adata: AnnData, name: str,
    *,
    per: Optional[str] = None,
    must_exist: bool = False,
) -> None:
    '''
    Delete a per-observation-per-variable matrix in ``adata`` by its ``name``.

    If ``per`` is specified, it must be one of ``obs`` (cells) or ``var`` (genes). This will only
    delete the cached layout-specific data. If ``per`` is not specified, both the data and any
    cached layout-specific data will be deleted.

    If ``must_exist``, will ``raise`` if the data does not currently exist.

    .. note::

        You can't delete the layer of the ``X`` data.

    .. todo::

        The function :py:func:`del_data_layer` is mainly required due to the need to delete cached
        layout-specific data. A better way would be to intercept ``del`` of the ``AnnData``
        ``layers`` field, but we do not own it.
    '''
    assert name.startswith('__')
    assert name != adata.uns['x_layer']

    if per is not None:
        per_name = '__per_%s__%s' % (per, name)
        if per_name in adata.layers:
            del adata.layers[per_name]
        elif must_exist:
            assert per_name in adata.layers
        return

    if must_exist:
        assert name in adata.layers
    for prefix in ['', '__per_obs__', '__per_var__']:
        prefixed_name = prefix + name
        if prefixed_name in adata.layers:
            del adata.layers[prefixed_name]


@utd.expand_doc()
def get_log_data(
    adata: AnnData,
    *,
    of_layer: Optional[str] = None,
    to_layer: Optional[str] = None,
    base: Optional[float] = None,
    normalization: float = 1,
    inplace: bool = False,
) -> utc.Matrix:
    '''
    Return a matrix with the log of some data layer.

    If ``of_layer`` is specified, this specific data is used. Otherwise, the data layer of ``X`` is
    used.

    If ``inplace``, the data will be stored in ``to_layer`` (by default, this is the name of the
    source layer with a ``log_`` prefix).

    The natural logarithm is used by default. Otherwise, the ``base`` is used.

    The ``normalization`` (default: {normalization}) is added to the count before the log is
    applied, to handle the common case of sparse data.

    .. note::

        The result is always a dense matrix, as even for sparse data, the log is rarely zero.
    '''
    if of_layer is None:
        of_layer = get_annotation_of_data(adata, 'x_layer')
        assert of_layer is not None

    if to_layer is None:
        to_layer = 'log_' + of_layer

    @timed.call()
    def compute_log_data() -> np.ndarray:
        assert of_layer is not None
        matrix = \
            get_data_layer(adata, of_layer, inplace=inplace)
        return utc.log_matrix(matrix, base=base, normalization=normalization)

    return get_data_layer(adata, to_layer, compute_log_data, inplace=inplace)


def get_annotation_of_obs(
    adata: AnnData,
    name: str,
    compute: Optional[Callable[[], utc.Vector]] = None,
    *,
    inplace: bool = False,
) -> np.ndarray:
    '''
    Lookup per-observation (cell) annotation in ``adata``.

    If the annotation does not exist, ``compute`` it. If no ``compute`` function was given,
    ``raise``.

    If ``inplace``, store the annotation in ``adata``.
    '''
    data = adata.obs.get(name)
    if data is not None:
        return data

    if compute is None:
        raise RuntimeError('unavailable observation annotation: ' + name)

    data = compute()
    assert data is not None
    assert data.shape == (adata.n_obs,)

    if inplace:
        utc.freeze(data)
        adata.obs[name] = data

    return data


def get_sum_per_obs(
    adata: AnnData,
    layer: str,
    *,
    inplace: bool = False,
    intermediate: bool = False,
) -> np.ndarray:
    '''
    Return the sum of the values per observation (cell) of some ``layer``.

    Use the existing ``sum_<layer>`` annotation if it exists.

    If ``inplace``, store the annotation in ``adata``.

    If ``intermediate``, store the layout-specific data used for the computation.
    '''
    @timed.call()
    def compute_sum_per_obs() -> np.ndarray:
        matrix = get_data_layer(adata, layer, per='obs', inplace=intermediate)
        return utc.sum_matrix(matrix, axis=1)

    return get_annotation_of_obs(adata, 'sum_' + layer, compute_sum_per_obs, inplace=inplace)


def get_nnz_per_obs(
    adata: AnnData,
    layer: str,
    *,
    inplace: bool = False,
    intermediate: bool = False
) -> np.ndarray:
    '''
    Return the number of non-zero values per observation (cell) of some ``layer``.

    Use the existing ``nnz_<layer>`` annotation if it exists.

    If ``inplace``, store the annotation in ``adata``.

    If ``intermediate``, store the layout-specific data used for the computation.
    '''
    @timed.call()
    def compute_nnz_per_obs() -> np.ndarray:
        matrix = get_data_layer(adata, layer, per='obs', inplace=intermediate)
        return utc.nnz_matrix(matrix, axis=1)

    return get_annotation_of_obs(adata, 'nnz_' + layer, compute_nnz_per_obs, inplace=inplace)


def get_max_per_obs(
    adata: AnnData,
    layer: str,
    *,
    inplace: bool = False,
    intermediate: bool = False
) -> np.ndarray:
    '''
    Return the maximal value per observation (cell) of some ``layer``.

    Use the existing ``max_<layer>`` annotation if it exists.

    If ``inplace``, store the annotation in ``adata``.

    If ``intermediate``, store the layout-specific data used for the computation.
    '''
    @timed.call()
    def compute_max_per_obs() -> np.ndarray:
        matrix = get_data_layer(adata, layer, per='obs', inplace=intermediate)
        return utc.max_matrix(matrix, axis=1)

    return get_annotation_of_obs(adata, 'max_' + layer, compute_max_per_obs, inplace=inplace)


def get_annotation_of_var(
    adata: AnnData,
    name: str,
    compute: Optional[Callable[[], utc.Vector]] = None,
    *,
    inplace: bool = False
) -> np.ndarray:
    '''
    Lookup per-variable (gene) annotation in ``adata``.

    If the annotation does not exist, ``compute`` it. If no ``compute`` function was given,
    ``raise``.

    If ``inplace``, store the annotation in ``adata``.
    '''
    data = adata.var.get(name)
    if data is not None:
        return data

    if compute is None:
        raise RuntimeError('unavailable variable annotation: ' + name)

    data = compute()
    assert data is not None
    assert data.shape == (adata.n_vars,)

    if inplace:
        utc.freeze(data)
        adata.var[name] = data

    return data


def get_sum_per_var(
    adata: AnnData,
    layer: str,
    *,
    inplace: bool = False,
    intermediate: bool = False
) -> np.ndarray:
    '''
    Return the sum of the values per variable (gene) of some ``layer``.

    Use the existing ``sum_<layer>`` annotation if it exists.

    If ``inplace``, store the annotation in ``adata``.
    '''
    @timed.call()
    def compute_sum_per_var() -> np.ndarray:
        matrix = get_data_layer(adata, layer, per='var', inplace=intermediate)
        return utc.sum_matrix(matrix, axis=0)

    return get_annotation_of_var(adata, 'sum_' + layer, compute_sum_per_var, inplace=inplace)


def get_nnz_per_var(
    adata: AnnData,
    layer: str,
    *,
    inplace: bool = False,
    intermediate: bool = False
) -> np.ndarray:
    '''
    Return the number of non-zero values per variable (gene) of some ``layer``.

    Use the existing ``nnz_<layer>`` annotation if it exists.

    If ``inplace``, store the annotation in ``adata``.

    If ``intermediate``, store the layout-specific data used for the computation.
    '''
    @timed.call()
    def compute_nnz_per_var() -> np.ndarray:
        matrix = get_data_layer(adata, layer, per='var', inplace=intermediate)
        return utc.nnz_matrix(matrix, axis=0)

    return get_annotation_of_var(adata, 'nnz_' + layer, compute_nnz_per_var, inplace=inplace)


def get_max_per_var(
    adata: AnnData,
    layer: str,
    *,
    inplace: bool = False,
    intermediate: bool = False
) -> np.ndarray:
    '''
    Return the maximal value per variable (gene) of some ``layer``.

    Use the existing ``max_<layer>`` annotation if it exists.

    If ``inplace``, store the annotation in ``adata``.

    If ``intermediate``, store the layout-specific data used for the computation.
    '''
    @timed.call()
    def compute_max_per_var() -> np.ndarray:
        matrix = get_data_layer(adata, layer, per='var', inplace=intermediate)
        return utc.max_matrix(matrix, axis=0)

    return get_annotation_of_var(adata, 'max_' + layer, compute_max_per_var, inplace=inplace)


@timed.call()
def get_annotation_of_data(
    adata: AnnData,
    name: str,
    compute: Optional[Callable[[], utc.Vector]] = None,
    *,
    inplace: bool = False
) -> Any:
    '''
    Lookup unstructured (data) annotation in ``adata``.

    If the annotation does not exist, ``compute`` it. If no ``compute`` function was given,
    ``raise``.

    If ``inplace``, store the annotation in ``adata``.
    '''
    data = adata.uns.get(name)
    if data is not None:
        return data

    if compute is None:
        raise RuntimeError('unavailable unstructured annotation: ' + name)

    data = compute()
    assert data is not None

    if inplace:
        adata.uns[name] = data

    return data
