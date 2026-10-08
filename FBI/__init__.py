"""
Environment that has to be set before numpy or matplotlib are imported anywhere
"""

import ctypes
import os
import sys

# One BLAS thread per process. Workers inherit an already initialised BLAS when they are
# forked, so this can't be set later, and `cores` workers each running a full thread pool
# would oversubscribe the machine. Raise FBI_BLAS_THREADS if running with cores=1.
_blas_threads = os.environ.get('FBI_BLAS_THREADS', '1')
for _var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_var, _blas_threads)

# If numpy is imported before somehow
if sys.platform == 'darwin' and os.environ['VECLIB_MAXIMUM_THREADS'] == '1':
    try:
        ctypes.CDLL('/System/Library/Frameworks/Accelerate.framework/Accelerate').BLASSetThreading(1)  # Single threaded
    except (OSError, AttributeError):
        pass

os.environ.setdefault('MPLBACKEND', 'Agg')

del ctypes, os, sys, _var, _blas_threads
