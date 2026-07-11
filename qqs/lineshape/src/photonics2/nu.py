"""DEPRECATED/DORMANT — outer product helper.

Status: dormant. Not maintained, not imported by any active code path.
The module-level print statements were removed 2026-07-11 (they printed
test data on import, which violates Python import conventions).

Originally written by original repo author (87930).
"""
import numpy as np

def direct_mul(a, b):
    la = len(a)
    lb = len(b)
    return np.matmul(a.reshape(la, 1), b.reshape(1, lb))