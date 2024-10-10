# solver

from oncophylo.tl.solver._condor import ConDoR
from oncophylo.tl.solver._scite import SCITE
from oncophylo.tl.solver._infscite import infSCITE
from oncophylo.tl.solver._scorchard import scOrchard
from oncophylo.tl.solver._sphyr import SPhyR
from oncophylo.tl.solver._huntress import HUNTRESS
from oncophylo.tl.solver._ltorchard import ltOrchard
from oncophylo.tl.solver._scarlet import scarlet

__all__ = (ConDoR, HUNTRESS, infSCITE, ltOrchard, scarlet, SCITE, scOrchard, SPhyR)