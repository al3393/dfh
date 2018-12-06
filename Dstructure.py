'''Defines type D structures.'''

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

class DGenerator(Generator):
    ''' Represents a generator of Type D structure. 
    Distinguished by (python) identity. '''
    def __init__(self, parent, idem):
        ''' Every generator must have an idempotent. '''
        Generator.__init__(self,parent)
        self.idem = idem # ask akram # left idempotent CB and fill in 
    
    def toSimpleDGenerator(self, name):
        ''' Convert to a SimpleDGenerator with the given name.
        All fields are preserved, except `name` which is overwritten,
        and hash_val which is removed, if present.'''
        new_obj = SimpleDGenerator(self.parent, self.idem, name)
        new_obj.__dict__.update(self.__dict__)
        new_obj.name = name # to overwrite original name
        if hasattr(new_obj, '_hash_val'):
            del new_obj._hash_val # reset hash value
        return new_obj

class SimpleDGenerator(DGenerator, NamedObject):
    def __init__(self, parent, idem, name): # idem? ask akram, left
        ''' Specifies name in addition.'''
        DGenerator.__init__(self, parent, idem)
        NamedObject.__init__(self, name)
    
class DStructure(FreeModule): # CT(Ti)
    '''Represents a type D structure. Note delta() returns an element in the 
    tensor module Tensor((A,D))'''
    def __init__(self, ring, algebra, side):
        '''Specifies the algebra and the side of D action'''
        FreeModule.__init__(self, ring, algebra, side)
        assert isinstance(algebra, DGAlgebra)
        self.algebra = algebra
        self.side = side # left
        # Construct A tensor M ( A is A(dLT), M is CT(T) for our case)
        # tensor product
        
        # tuple remove
        def _mul_A_AtensorM(AGenMGen, ACoeff):
            '''To be used as rmultiply() in AtensorM. Multiply ACoeff with
            AGen. `AGenMGen` is a tuple of (AGen, MGen),'''
            AGen, MGen = AGenMGen
            return (ACoeff * AGen) * MGen

        def _diff_AtensorM(AGenMGen):
            '''To be used as diff() in AtensorM. AGenMGen` is a tuple of (AGen, MGen)'''
            AGen, MGen = AGenMGen
            return (AGen.diff() * MGen) + (AGen * MGen.delta())
        
        self.AtensorM.rmultiply = _mul_A_AtensorM # store the functions
        self.AtensorM.diff = _diff_AtensorM
      
    def delta(self, generator):
        ''' Returns delta^1 of the generator.'''
        raise NotImplementedError("Differential not implemented.")
    
    def rmultiply(self, MGen, AGen):
        ''' Multiply a generator of DStructure with an algebra generator
        means forming a tensor'''
        return 1 * TensorGenerator((AGen, MGen),self.AtensorM )
        
class SimpleDStructure(DStructure):
    ''' Represents a type D Structure with a finite number of generators, and 
    explicitly stored generating set and delta option'''
    
    def __init__(self, ring, algebra, side = ACTION_LEFT):
        ''' Initialize an empty type D Structure.'''
        assert side == ACTION_LEFT, "Right action not implemented."
        DStructure.__init__(self, ring, algebra, side)
        self.generators = set()
        self.delta_map = dict()
    
    def __len__(self):
        return len(self.generators)
    
    def delta(self, generator):
        return self.delta_map[generator]
    
    def getGenerators(self):
        return list(self.generators)
    
    def addGenerator(self, generator):
        '''Add a generator. No effect if the generator already exists. '''
        assert generator.parent == self
        assert isinstance(generator, DGenerator)
        self.generators.add(generator)
        if not generator in self.delta_map:
            self.delta_map[generator] = E0 # CB
    
    def addDelta(self, gen_from, gen_to, alg_coeff, ring_coeff):
        ''' Add ring_coeff(F2) * alg_coeff(aij) * gen_to(x) to the 
        delta of gen_From. Both arguments should generators'''
        assert gen_from.parent == self and gen_to.parent == self
        if alg_coeff is None:
            alg_coeff = gen_to.idem.toAlgElt(self.algebra)
        # CB and implement idempotent in A(dLT)
        assert alg_coeff.getLeftIdem() == gen_from.idem # CB and ask Akram
        assert alg_coeff.getRightIdem() == gen_to.idem # CB and ask Akdram - is it not del? is it DL?
        self.delta_map[gen_from] += (alg_coeff * gen_to) * ring_coeff
        
    def reindex(self):
        """Replace the generators by simple generators indexed by integers."""
        gen_list = list(self.generators)
        new_gen_list = []
        translate_dict = dict()
        for i in range(len(gen_list)):
            new_gen = gen_list[i].toSimpleDGenerator("g%d"%(i+1))
            new_gen_list.append(new_gen)
            translate_dict[gen_list[i]] = new_gen
        self.generators = set(new_gen_list)
        new_delta = dict()
        for k, v in self.delta_map.items():
            new_v = E0
            for (AGen, MGen), coeff in v.items(): # cb and check (AGen, MGen)
                new_v += (AGen * translate_dict[MGen]) * coeff
            new_delta[translate_dict[k]] = new_v
        self.delta_map = new_delta
        if hasattr(self, "grading"):
            new_grading = dict()
            for gen, gr in self.grading.items():
                if gen in translate_dict: # gen is still in dstr
                    new_grading[translate_dict[gen]] = gr
            self.grading = new_grading
    
    def deltaCoeff(self, gen_from, gen_to):
        ''' Return the coefficient ( as algebra element) of gen_to in delta of gen_from'''
        if self.delta_map[gen_from] == 0:
            return E0
        else:
            return self.delta_map[gen_from].fixLast(gen_to)
    
    def testDelta(self):
        ''' Verify d^2 for this structure. '''
        for gen in self.generators:
            if gen.delta().diff() != 0 or gen.delta().diff() != E0 :
                # print the terms in d^2 for one gerator
                print("{0} ==> ".format(gen))
                for k, v in gen.delta().diff().items():
                    print("{0} * {1}".format(v,k))
                return False
        return True
    
    def __str__(self):
        result = "Type D Structure. \n"
        for k,v in self.delta_map.items():
            result += "d(%s) = %s\n" % (k, v)
        return result
    
    def simplify(self, cancellation_constraint = None):
        """Simplify a type D structure using cancellation lemma."""
        # Simplification is best done in terms of coefficients
        # Build dictionary of coefficients
        arrows = dict()
        for gen in self.generators:
            arrows[gen] = dict()
        for gen in self.generators:
            for (AGen, MGen), coeff in self.delta_map[gen].items():
                if MGen not in arrows[gen]:
                    arrows[gen][MGen] = E0
                arrows[gen][MGen] += AGen * coeff

        arrows = simplifyComplex(
            arrows, E0,
            cancellation_constraint = cancellation_constraint)

        # Now rebuild the type D structure
        self.generators = set()
        self.delta_map = dict()
        for x in arrows:
            self.generators.add(x)
            self.delta_map[x] = E0
            for y, coeff in arrows[x].items():
                self.delta_map[x] += coeff * y

        # This is a good place to simplify gradings #CB and remove
        if hasattr(self, "gr_set"):
            new_gr_set = self.gr_set.simplifiedSet()
            for gen in self.generators:
                self.grading[gen] = self.gr_set.simplifiedElt(self.grading[gen])
            self.gr_set = new_gr_set
    
    def checkGrading(self):
        ''' Check grading is consistent with the type D operations.'''
        for x in self.generators:
            for (a,y), coeff in x.delta().items():
                gr_x = self.grading[x] # cb and implement grading
                gr_y = self.grading[y]
                assert gr_x[0] - 1 == gr_y[0] + a.getGrading()[0] # maslov is maintained
                assert gr_x[1] == gr_y[0] + a.getGradig()[1] # alexander is maintained
                           #CB akram
    def compareDStructures(self, other):
        """Compare two type D structures, print out any differences."""
        # Some basic tests:
        if len(self) != len(other):
            print("Different number of generators.")
            return False
        if self.algebra != other.algebra:
            print("Different algebra action.")
            return False
        gen_map = dict()
        for gen1 in self.generators:
            for gen2 in other.generators:
                if gen1.idem == gen2.idem:
                    gen_map[gen1] = gen2
                    break
        for gen1 in self.generators:
            for gen2 in self.generators:
                coeff1 = self.deltaCoeff(gen1, gen2)
                coeff2 = other.deltaCoeff(gen_map[gen1], gen_map[gen2])
                if coeff1 != coeff2:
                    print("Different coefficient at %s->%s" % (gen1, gen2))
                    print("%s vs %s" % (coeff1, coeff2))
                    return False
        return True