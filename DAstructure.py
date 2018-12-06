'''Defines type DA Structures.'''

from algebra import DGAlgebra, FreeModule, Generator, Tensor, \
    TensorGenerator, TensorStarGenerator, AGenerator
from algebra import ChainComplex, E0, TensorDGAlgebra
from dstructure import DGenerator, SimpleDStructure
from grading import GeneralGradingSet, GeneralGradingSetElement, \
    SimpleDbGradingSetElement
from hdiagram import getIdentityDiagram
from pmc import StrandDiagram
from utility import MorObject, NamedObject
from utility import sumColumns
from utility import ACTION_LEFT, ACTION_RIGHT, F2

class DAGenerator(Generator):
    ''' Represents a generator of type DA Structure.'''
    def __init__(self, parent, idem1, idem2):
        ''' Every generator has two idempotents. 
        idem1: left type D idempotent on the left.
        idem2: right type A idempotent on the right.''' # cb and check
        Generator.__init__(self,parent)
        self.idem1, self.idem2 = idem1, idem2
    
class SimpleDAGenerator(DAGenerator, NamedObject):
    ''' Represents a generator of type DA structure, distinguished by name. '''
    def __init__(self, parent, idem1, idem2, name):
        ''' Specifies name in addition.'''
        DAGenerator.__init__(self, parent, idem1, idem2)
        NamedObject.__init__(self,name)
    
    def __str__(self):
        return "%s:%s,%s" % (self.name, self.idem1, self.idem2)
    
    def __repr__(self):
        return str(self)

class DTensorAenerator(DAGenerator, tuple): # cb and fill in
    pass

class DATensorDGenerator(DGenerator, tuple): # cb and fill in
    ''' Generator of a type D structure formed by tensoring 
    a type DA structure and a type D structure. ''' # akram
    
    def __new__(cls, parent, gen_left, gen_right):
        ''' gen_left is a generator of a type DA structure.
        gen_right is a generator of type D structure.'''
        return tuple.__new__(cls, (gen_left, gen_right))
    
    def __init__(self, parent, gen_left, gen_right):
        DGenerator.__init__(self, parent, gen_left.idem) # missing come back and fill in

class ATensorDGenerator(AGenerator, tuple):
    ''' Generator of a type A structure formed by tensoring 
    a type A structure and a type DA structure. '''
    pass

class DATensorDAGenerator(DAGenerator, tuple):
    pass
    
class DAStructure(FreeModule):
    ''' Represents a type DA structure. delta() takes a generator and a sequence
    of algebra generators (as a generator of the tensor algebra) and returns
    an element in the tensor module tensor ((A,M)), where A is algebra (D-side) algebra'''
    
    def __init__(self, ring, algebra1, algebra2, side1, side2):
        '''Specifies the algebras and sides of the type DA action'''
        FreeModule.__init__(self, ring)
        assert isinstance(algebra1, DGAlgebra)
        assert isinstance(algebra2, DGAlgebra)
        self.algebra1 = algebra1
        self.side = side1
        self.algebra2 = algebra2
        self.side2 = side2
        
        