# solver

from oncophylo.tl.solver._condor import ConDoR
from oncophylo.tl.solver._compass import COMPASS
from oncophylo.tl.solver._scite import SCITE
from oncophylo.tl.solver._infscite import infSCITE
from oncophylo.tl.solver._scaorchard import scaOrchard
from oncophylo.tl.solver._sphyr import SPhyR
from oncophylo.tl.solver._huntress import HUNTRESS
from oncophylo.tl.solver._ltorchard import ltOrchard
from oncophylo.tl.solver._scarlet import scarlet
from oncophylo.tl.solver._grmt import GRMT
from oncophylo.tl.solver._bitsc2 import BitSC2

__all__ = (BitSC2, ConDoR, COMPASS, GRMT, HUNTRESS, infSCITE, ltOrchard, scarlet, SCITE, scaOrchard, SPhyR)