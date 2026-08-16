"""
Environment that has to be set before numpy or matplotlib are imported anywhere
"""

import os

# One BLAS thread per process. Workers inherit an already initialised BLAS when they are
# forked, so this can't be set later, and `cores` workers each running a full thread pool
# would oversubscribe the machine. Raise FBI_BLAS_THREADS if running with cores=1.
_blas_threads = os.environ.get('FBI_BLAS_THREADS', '1')
for _var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_var, _blas_threads)

# pydarn imports matplotlib, and the macOS backend starts Objective-C runtime state that
# does not survive a fork
os.environ.setdefault('MPLBACKEND', 'Agg')

del os, _var, _blas_threads
