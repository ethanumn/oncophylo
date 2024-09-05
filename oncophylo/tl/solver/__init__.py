# solver

from oncophylo.tl.solver._condor import ConDoR
from oncophylo.tl.solver._scite import SCITE
from oncophylo.tl.solver._infscite import infSCITE
from oncophylo.tl.solver._scorchard import scOrchard
from oncophylo.tl.solver._sphyr import SPhyR

__all__ = (ConDoR, infSCITE, SCITE, scOrchard, SPhyR)