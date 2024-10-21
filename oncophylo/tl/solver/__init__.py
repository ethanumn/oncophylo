# solver

from oncophylo.tl.solver._condor import ConDoR
from oncophylo.tl.solver._scite import SCITE
from oncophylo.tl.solver._infscite import infSCITE
from oncophylo.tl.solver._scorchard import scOrchard
from oncophylo.tl.solver._sphyr import SPhyR
from oncophylo.tl.solver._huntress import HUNTRESS
from oncophylo.tl.solver._ltorchard import ltOrchard
from oncophylo.tl.solver._scarlet import scarlet
from oncophylo.tl.solver._grmt import GRMT

__all__ = (ConDoR, GRMT, HUNTRESS, infSCITE, ltOrchard, scarlet, SCITE, scOrchard, SPhyR)