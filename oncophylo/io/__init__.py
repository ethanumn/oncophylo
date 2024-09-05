"""I/O"""

from oncophylo.io._genotype import load
from oncophylo.io._tree import load_dot


__all__ = (load, load_dot)