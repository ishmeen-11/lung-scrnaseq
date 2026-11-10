"""lungsc: preprocessing for the Prince lab 3-species lung scRNA-seq data."""
__version__ = "0.1.0"

# Compatibility shim: newer SciPy dropped .getnnz() on sparse *arrays*, but
# scanpy's Scrublet/neighbors code still calls it. Restore it as an alias.
from scipy import sparse as _sp

for _cls in (_sp.csr_array, _sp.csc_array, _sp.coo_array):
    if not hasattr(_cls, "getnnz") and hasattr(_cls, "_getnnz"):
        _cls.getnnz = lambda self, axis=None: self._getnnz(axis=axis)
