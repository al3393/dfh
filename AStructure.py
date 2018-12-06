'''Defines type A structures.'''

from fractions import Fraction
from algebra import DGAlgebra, FreeModule, Generator, SimpleChainComplex, \
    Tensor, TensorGenerator
from algebra import simplifyComplex
from algebra import E0
from grading import GeneralGradingSet, GeneralGradingSetElement
#from hdiagram import getZeroFrameDiagram, getInfFrameDiagram, getPlatDiagram
from tangle_new import Idempotent, Strands, StrandDiagram
from pmc import connectSumPMC, splitPMC, linearPMC
from utility import MorObject, NamedObject
from utility import memorize
from utility import ACTION_LEFT, DEFAULT_GRADING, F2, SMALL_GRADING

class SimpleAStructure(AStructure):
    ''' Represents a type D Structure with a finite number of generators, and 
    explicitly stored generating set and delta option'''
    
    def __init__(self, ring, algebra, side = ACTION RIGHT):
        ''' Initialize an empty type D Structure.'''
        assert side == ACTION_LEFT, "Right action not implemented."
        DStructure.__init__(self, ring, algebra, side)
        self.generators = set()
        self.delta_map = dict()
    
    def addm_2(self, gen_from, gen_to, alg_coeff, ring_coeff):
        ''' Add ring_coeff(F2) * alg_coeff(aij) * gen_to(x) to the 
        delta of gen_From. Both arguments should generators'''
        assert gen_from.parent == self and gen_to.parent == self
        if alg_coeff is None:
            alg_coeff = gen_to.idem.toAlgElt(self.algebra)
        # CB and implement idempotent in A(dLT)
        assert alg_coeff.getLeftIdem() == gen_from.idem # CB and ask Akram
        assert alg_coeff.getRightIdem() == gen_to.idem # CB and ask Akdram - is it not del? is it DL?
        self.delta_map[gen_from] += (alg_coeff * gen_to) * ring_coeff