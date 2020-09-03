'''
Use multiple CPUs.

Due to the notorious GIL, using multiple Python threads is useless. This leaves us
with two options for using multiple CPUs:

* Use multiple threads in the internal C++ implementation of some Python functions;
  this is done by both numpy and the C++ extension functions provided here, and works for reasonably
  small sized work, such as sorting each of the rows of a large matrix.

* Use Python multi-processing. This is costly and only works for large sized work, such as
  computing metacells for different piles.

Each of these two approaches works tolerably well on its own, even though both are sub-optimal. The
problem starts when we want to combine them. Consider a server with 50 processors. Invoking
``corrcoef`` on a large matrix will use them all. This is great if one computes metacells for a
single pile. Suppose, however, you want to compute metacells for 50 piles, and do so using
multi-processing. Each and every of the 50 sub-processes will invoke ``corcoeff`` which will spawn
50 internal threads, resulting in the operating system seeing 2500 processes competing for the same
50 hardware processors. This "does not end well".

The proper solution would be to get rid of the GIL, not use multi-processing, only use multiple
threads in a pool, and adopt a universal mechanism for running tasks on the threads pool. This will
ensure we have at most 50 threads, used even if parallel tasks spawn nested parallel tasks, and all
memory would be shared so there would be zero-cost in accessing large data structures from different
threads.

You would expect that, a decade after multi-core CPUs became available, this would have been solved
"out of the box", by the current frameworks agreeing to cooperate with each other. However, somehow
this isn't seen as important by the people maintaining these frameworks; in fact, most of them don't
properly handle nested parallelism within their own framework, never mind playing well with others.

So in practice, while new languages (such as Julia and Rust) provide this out of the box, old
languages (such as Python and C++) are stuck in a swamp. In our case, numpy uses some underlying
parallel threads framework, our own extensions use OpenMP parallel threads, and we are forced to use
the Python-multi-processing framework itself oname top of both, where each of these frameworks is blind
to the rest.

As a crude band-aid, there are ways to force both whatever-numpy-uses and OpenMP to use a specific
number of threads. So, when we use multi-processing, we limit each sub-process to use less internal
threads, such that the total will be at most 50. This is even less optimal, but at least it doesn't
bring the server to its knees trying to deal with a total load of 2500 processes.

.. note::

    Even if only using internal threads, it seems that the total CPU time grows significantly when
    using multiple threads. The total elapsed time is reduced, but since the total is increased,
    the gain isn't nearly as much as one would have hoped.

.. todo::

    Look into using TBB instead of OpenMP for the inner threads for better efficiency.
'''
import ctypes
import logging
import os
import sys
from multiprocessing import Pool, Value
from threading import current_thread
from typing import Any, Callable, Iterable, Optional

from threadpoolctl import threadpool_limits  # type: ignore

__all__ = [
    'set_cpus_count',
    'get_cpus_count',
    'parallel_map',
]


LOG = logging.getLogger(__name__)


CPUS_COUNT = 0

MAIN_PROCESS_PID = os.getpid()
IS_MAIN_PROCESS: Optional[bool] = True

PROCESS_INDEX = 0
current_thread().name = '#0'

PROCESSES_COUNT = 0
NEXT_PROCESS_INDEX = Value(ctypes.c_int32, lock=True)
PARALLEL_FUNCTION: Optional[Callable[[int], Any]] = None


def set_cpus_count(cpus: int) -> None:
    '''
    Set the (maximal) number of CPUs to use in parallel.

    By default, we use all the available hardware threads. Override this by setting the
    ``METACELLS_CPUS_COUNT`` environment variable or by invoking this function from the main thread.

    A value of ``-1`` uses half of the available threads to combat hyper-threading, ``0`` uses all
    the available threads, and otherwise the value is the number of threads to use.
    '''
    assert IS_MAIN_PROCESS

    if cpus == -1:
        cpus = (os.cpu_count() or 1) // 2
    elif cpus == 0:
        cpus = os.cpu_count() or 1

    assert cpus > 0

    global CPUS_COUNT
    CPUS_COUNT = cpus

    threadpool_limits(limits=CPUS_COUNT)


if not 'sphinx' in sys.argv[0]:
    set_cpus_count(int(os.environ.get('METACELLS_CPUS_COUNT',
                                      str(os.cpu_count()))))


def get_cpus_count() -> int:
    '''
    Return the number of CPUs we are allowed to use.
    '''
    assert CPUS_COUNT > 0
    return CPUS_COUNT


def parallel_map(
    function: Callable[[int], Any],
    invocations: int,
) -> Iterable[Any]:
    '''
    Execute ``function``, in parallel, ``invocations`` times. Each invocation is given the
    invocation's index as its single argument.

    For our simple pipelines, only the main process is allowed to execute functions in parallel
    processes, that is, we do not support nested ``parallel_map`` calls.

    If the current :py:func:`get_cpus_count` is one, just runs the function serially. Otherwise,
    fork new processes to execute the function invocations (using ``multiprocessing.Pool``). The
    downside is that this is slow. The upside is that each of these processes starts with a shared
    memory copy(-on-write) of the full Python state, that is, all the inputs for the function are
    available "for free".
    '''
    assert function.__is_timed__  # type: ignore

    global IS_MAIN_PROCESS
    assert IS_MAIN_PROCESS

    global PROCESSES_COUNT
    PROCESSES_COUNT = min(CPUS_COUNT, invocations)
    if PROCESSES_COUNT == 1:
        return [function(index) for index in range(invocations)]

    NEXT_PROCESS_INDEX.value = 0

    global PARALLEL_FUNCTION
    assert PARALLEL_FUNCTION is None

    PARALLEL_FUNCTION = function
    IS_MAIN_PROCESS = None
    try:
        with Pool(PROCESSES_COUNT) as pool:
            return pool.map(_invocation, range(invocations))
    finally:
        IS_MAIN_PROCESS = True
        PARALLEL_FUNCTION = None


def _invocation(index: int) -> Any:
    global IS_MAIN_PROCESS
    if IS_MAIN_PROCESS is None:
        IS_MAIN_PROCESS = os.getpid() == MAIN_PROCESS_PID
        assert not IS_MAIN_PROCESS

        global PROCESS_INDEX
        global NEXT_PROCESS_INDEX
        with NEXT_PROCESS_INDEX:  # type: ignore
            PROCESS_INDEX = NEXT_PROCESS_INDEX.value
            NEXT_PROCESS_INDEX.value += 1

        current_thread().name = '#%s' % PROCESS_INDEX

        global CPUS_COUNT
        start_cpu_index = \
            round(CPUS_COUNT * PROCESS_INDEX / PROCESSES_COUNT)
        stop_cpu_index = \
            round(CPUS_COUNT * (PROCESS_INDEX + 1) / PROCESSES_COUNT)
        CPUS_COUNT = stop_cpu_index - start_cpu_index

        assert CPUS_COUNT > 0
        threadpool_limits(limits=CPUS_COUNT)

    assert PARALLEL_FUNCTION is not None
    return PARALLEL_FUNCTION(index)