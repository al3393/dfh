#from numbers import Number
#import heapq
#from numbers import Number
from utility import NamedObject, SummableDict,safeMultiply
#from utility import fracToInt, memorize, memorizeHash, safeMultiply
#from utility import F2
from numbers import Number
import heapq

class FreeModule:
    '''Represents a free module over some ring.'''
    
    def __init__(self,ring):
        '''Specifies the ring. Ring should be a subclass of Ring or a python number'''
        self.ring = ring
    
    def getGenerators(self):
        '''returns a list of all generators. Need not be implemented for every module.'''       
        raise NotImplementedError("Method getGenerators not implemented.")

class Generator:
    '''Represents a generator of a free module.
    Implement __eq__, __ne__, and __hash__.'''
    
    def __init__(self,parent):
        '''every Generator needs the info of its "parent" module'''
        self.parent = parent
    
    def elt(self, coeff = 1):
        '''Returns the element coeff * self'''
        return self.ELT_CLASS({self: coeff})  # returns the element(summable dict)
                                        # object with key self, value coeff
    def diff(self):
        '''returns the differential of this generator.
        It calls the diff() method of parent module.'''
        return self.parent.diff(self)
    
    def factor(self):
        ''' Finding all ways to factor this generator into a product of
        two generators. Calls factor() method of its parent module''' 
        return self.parent.factor(self)

    def delta(self):
        '''Returns the delta of this generator. 
        It calls the delta() method of parent module.'''      
        return self.parent.delta(self)
    
    def __mult__(self, other):
        '''multiples this generator by ``other`` on the left.
        That is, the right module action'. Usually returns an element
        rather than a generator'''   
        if isinstance(other,Number):
            return self.elt(other)
        elif hasattr(self.parent, "multiply"):
            return self.parent.multiply(self,other)
        else:
            return NotImplemented
    
    def __rmul__(self,other):
        '''multiplies this generator by ``other`` on the right. 
        That is, the left module ation
        Usually represents a left module action''' 
        if isinstance(other, Number):
            return self.elt(other)
        elif hasattr(self.parent, "rmultiply"):
            return self.parent.rmultiply(self, other)
        else:
            return NotImplemented

    def toSimpleGenerator(self, name): # ? Necessary? 
        '''Convert to a SimpleGenerator with the given name. All fields are
        preserved, except ``name`` which is overwritten, and _hash_val which is
        removed, if present.'''
        new_obj = SimpleGenerator(self.parent, name)
        new_obj.__dict__.update(self.__dict__)
        new_obj.name = name # to make sure original name is overwritten
        if hasattr(new_obj, '_hash_val'):
            del new_obj._hash_val # reset hash value
        return new_obj

class SimpleGenerator(Generator, NamedObject):
    '''Each Generator has a name. Distinguished by name.'''
    def __init__(self, parent, name):
        """Specifies name and parent module."""
        Generator.__init__(self, parent)
        NamedObject.__init__(self, name) # ?
    
class Element(SummableDict):
    '''Represents an element of a free module, as a dictionary from generators
    to coefficients. ex) a+3b will be written as {a:1,b:3}'''
    
    def __init__(self,data = None):
        if data is None:
            data = {}

        SummableDict.__init__(self,data)
        
        if self: # if the object is instantialized
            convert = self.getElt().parent.ring.convert #?
            for key,value in self.items():
                self[key] = convert(value)
    
    def __str__(self):
        #string method for Element Class
        if self == 0:
            return "0"
        terms = []
        for gen, coeff in self.items():
            if coeff == 1:
                terms.append(str(gen))
            else:
                terms.append(str(coeff)+"*"+str(gen))
        return "+".join(terms)
    
    def __repr__(self):
        return str(self)
    
    def __mul__(self, other):
        # First try multiplying each coefficient with other, using the function
        # in SummableDict.
        result = SummableDict.__mul__(self, other)
        if result != NotImplemented:
            return result

        # Now try to multiply each key by other on the right.
        result = E0
        for k, v in self.items():
            prod = safeMultiply(k, other)
            if prod is NotImplemented:
                return NotImplemented
            result += [term * (v * coeff) for term, coeff in prod.items()]
        return result
    
    def __rmul__(self, other):
        # First try multiplying each coefficient with other, using the function
        # in SummableDict.
        result = SummableDict.__rmul__(self, other)
        if result != NotImplemented:
            return result
        # Now try to multiply key by other on the left.
        result = E0
        for k, v in self.items():
            prod = safeMultiply(other, k)
            if prod is NotImplemented:
                return NotImplemented
            result += [term * (v * coeff) for term, coeff in prod.items()]
        return result
    
    def diff(self):
        ''' Returns the differential of this element.'''
        return sum([coeff * gen.diff() for gen, coeff in self.items()], E0)
  
# Name of the class for elements containing this generator
Generator.ELT_CLASS = Element
# Short-hand for empty element
E0 = Element()

class ChainComplex(FreeModule): # ask akram : why is this extension of free module
    '''Represents a general chain complex. '''
    def diff(self, gen):
        ''' Returns the differential of a generator.'''
        raise NotImplementedError("Differential not implemented.")

class SimpleChainComplex(ChainComplex): 
    '''Represents a chain complex with a finite number of generators, with 
    explicitly stored generating set and differential.Generating set is stored
    as a python set, and differential is stored as a dictionary mapping
    from generators to elements. Each generator must be a key in the dictionary
    (even if its differential is zero)''' 
    def __init__(self,ring):
        ''' Initialize an empty chain complex.'''
        ChainComplex.__init__(self,ring)
        self.generators = set()
        self.differential = dict()
    
    def __str__(self):
        result = "Chain complex. \n"
        for k,v in self.differential.items():
            result += "d(%s) = %s \n" %(k,v)
        return result

    def __repr__(self):
        return str(self)
    
    def __len__(self):
        return len(self.generators)
    
    def diff(self, generator):
        return self.differential[generator]
    
    def diffElt(self, elt):
        '''Return the differential of Element elt of this SimpleChainComplex'''
        # applying linearity to element
        answer = E0
        for x in elt.keys():
            answer += elt[x] * self.diff(x)
        return answer
    
    def getGenerators(self):
        return list(self.generators)
    
    def reindex(self): # cb and fill in
        '''Replace the generators by simple generators indexed by integers.
        The names of the new generators are 'g1', 'g2', etc.'''
        gen_list = list(self.generators)
        new_gen_list = []
        # Dictionary mapping original generators to new ones
        translate_dict = dict()
        for i in range(len(gen_list)):
            new_gen = gen_list[i].toSimpleGenerator("g%d"%(i+1))
            new_gen_list.append(new_gen)
            translate_dict[gen_list[i]] = new_gen
        self.generators = set(new_gen_list)
        new_diff = dict()
        for gen, dgen in self.differential.items():
            new_diff[translate_dict[gen]] = dgen.translateKey(translate_dict)
        self.differential = new_diff

        #translate Maslov grading CB and modify
        if hasattr(self, "m_grading"):
            new_grading = dict()
            for gen, gr in self.m_grading.items():
                if gen in translate_dict: # gen is still in chain complex
                    new_grading[translate_dict[gen]] = gr
            self.m_grading = new_grading
        
        #translate Alexander grading CB and modify
        if hasattr(self, "a_grading"):
            new_grading = dict()
            for gen, gr in self.a_grading.items():
                if gen in translate_dict: # gen is still in chain complex
                    new_grading[translate_dict[gen]] = gr
            self.a_grading = new_grading

    def addGenerator(self, generator):
        '''Add an generator. No effect if the generator already exists.'''
        assert generator.parent == self
        self.generators.add(generator)
        #if not self.differential.has_key(generator): # has key removed in python 3
        if generator not in self.differential:
            self.differential[generator] = E0    
    
    def addDifferential(self, gen_from, gen_to, coeff):
        ''' Add coeff * gen_to to the differential of gen_from.
        i.e diff(gen_from) = coeff * gen_to
        gen_from and gen_to must be generators of this complex
        '''
        assert gen_from.parent == self and gen_to.parent == self
        self.differential[gen_from] += coeff * gen_to
        
    #CB ask akram
    def simplify(self, find_homology_basis = False, cancellation_constraint = None): 
        # Build dictionary of coefficients
        arrows = dict()
        for x in self.generators:
            arrows[x] = dict()
        for x in self.generators():
            for y, coeff in self.differential[x].items():
                arrows[x][y] = coeff
        arrows =  simplifyComplex(arrows, find_homology_basis = True) # cb and change
        # CB and change

    def checkDifferential(self):
        '''Check the relations d^2 for differentials.'''
        for gen in self.generators:
            assert gen.diff().diff() == 0
    
    def checkGrading(self): # cb and fix grading[x] - where is the grading method?
        ''' Check if the grading is consistent with differentials
        i.e., Alexander is preserved, and Maslov down by one '''
        for x in self.generators:
            for y, coeff in x.diff().items():
                assert self.m_grading[x] - 1 == self.m_grading[y] \
                        and self.a_grading[x] == self.a_grading[y]
    
    def getGradingInfo(self): # CB and fill in - remove if unncessary
        ''' Shows the distribution of gradings in an easy-to-read format.'''
        distr_by_maslov = dict()
        distr_by_alex = dict()

    def copy(self): # cb and fill in 
        ''' Return a copy of `self'''
        new_to_old = dict()
        old_to_new = dict()
        answer = SimpleChainComplex(self.ring)
        genlist = list(self.generators)
        for i in range(len(self.generators)):
            newgen = SimpleGenerator(answer,repr(i))
            answer.addGenerator(newgen)
            new_to_old[newgen] = genlist[i]
            old_to_new[genlist[i]] = newgen
        for x in genlist:
            for y in self.differential[x].keys():
                answer.addDifferential(old_to_new[x],old_to_new[y],self.differential[x][y])
        return answer
        
    def id(self): # cb and remove if unnecessary ask akram
        pass

class DGAlgebra(ChainComplex):
    '''Represents a general differential-graded algebra.'''
    def multiply(self, gen1, gen2):
        '''Returns the product of gen1 and gen 2, as an algebra element.'''
        raise NotImplementedError("Multiply not implemented.")
    
    def factormap(): # CB and remove if not necessary
        pass
        
class Tensor(FreeModule, tuple):
    '''Represents a free module whose generating set is the product of the generating set of two sides.''' 
        #two sides ?  ask akram
    #tuples
    def __init__(self,data):
        '''Specifies the left and right module'''
        for d in data[1:]:
            assert d.ring == data[0].ring # raise assertion if the base ring is not the same
    #Note that the tuple initialization is automatic 
        FreeModule.__init__(self,data[0].ring) # specifies the generator and the rings 

class TensorGenerator(Generator, tuple): #?
    '''Represents a generator of a free module that is a tensor product two more free modules. 
    Works as a tuple of the components'''
    
    def __new__(cls, data, parent = None):
        return tuple.__new__(cls, tuple(data))
    
    def __init__(self, data, parent = None):
        ''' If parent is None, a default is used.'''
        if parent is None:
            parent = Tensor(tuple([comp.parent for comp in data]))
        # Note tuple initialization is automatic
        Generator.__init__(self, parent)
    
    def __eq__(self,other):
        return tuple.__eq__(self,other) and self.parent == other.parent
    
    def __ne__(self, other):
        return not (self == other)
    
    def __hash__(self):
        return hash((tuple(self), "Tensor"))
    
    def __str__(self):
        return "*".join(str(comp) for comp in self)
    
    def getLeftIdem(self): # CB and check
        ''' Get the left idempotent. Only works if same function is implemented 
        in each component.'''
        return TensorIdempotent(tuple([comp.getLeftIdem() for comp in self]))
    
    def getRightIdem(self): # CB and check
        ''' Get the right idempotent. Only works if same function is implemented 
        in each part.'''
        
        return TensorIdempotent(tuple([comp.getRightIdem() for comp in self]))
    
class TensorElement(Element):
    '''Represents an element of the tensor product of two or more modules.
    CB and finish the support quickly collecting terms by one of the components. '''
    
    def __init__(self, data= None, parent = None):
        ''' If the keys are tuples, convert them to tensor generators. '''
        if data is None:
            data = {}
        data_processed ={}
        #a = (('a',3),('b',4),('z',5)) converted to dic is {'a': 3, 'b': 4, 'z': 5}
        for term, coeff in dict(data).items():
            if isinstance(term, TensorGenerator):
                data_processed[term] = coeff
            else:
                data_processed[TensorGenerator(term,parent)] = coeff
        Element.__init__(self, data_processed)

    def fixLast(self, gen, parent = None):
        ''' Collect terms with `gen` as the last factor. Return the coefficient as
        either Element or Tensor Element'''
        # Returning in element form 
        result = E0 
        for term, coeff in self.items():
            if term[-1] == gen: # if last term is gen
                if len(term) > 2:
                    result += coeff * TensorGenerator(term[:-1], parent) # omit last character
                else: # len(term) == 2
                    result += coeff * term[0]
        return result
    
    def invertible(self):
        ''' Tests whether this element is invertible. '''
        for term, coeff in self.items():
            for comp in term:
                if not (1*comp).invertible():
                    return False
        return True
    
    def inverse(self):
        ''' Returns the inverse of this element, if invertible. Undefined behavior
        if the element is not invertible'''
        return self
    
TensorGenerator.ELT_CLASS = TensorElement

# applying linearity to tensor product
def expandTensor(prod, parent = None):
    '''Produces the tensor element formed by the tensor product of either
    generators or elements.
    
    ``prod`` is a tuple of either Generator or Element, corresponding to the
    components of the tensor product.
    
    For example, ((1*A+1*B),C) expands into 1*(A,C)+1*(B,C), and
    ((1*A-1*B),(1*C-1*D)) expands into 1*(A,C)-1*(B,C)-1*(A,D)+1*(B,D).
    ASSOCIATIVITY/ LINEARITY ( A,B,C,D are elements)
    
    ``parent`` specifies the Tensor module. If it is set to None, the default
    (with no additional operations defined) will be used (during the
    initialization of TensorElement).'''
    
    assert isinstance(prod, tuple)
    num_part = len(prod)
    expanded = [(prod,1)]
    for i in range(num_part):
        if len(expanded) == 0:
            return E0
        if isinstance(expanded[0][0][i],Generator):
            continue
        expanded2 = []
        for subterm, coeff in expanded:
            for gen, coeff2 in subterm[i].items():
                expanded2.append(((subterm[0:i]+(gen,)+subterm[i+1:]),
                                 coeff*coeff2))
        expanded = expanded2
    
    if isinstance(parent, TensorStar): # if parent is Tensor star.
    # i.e a free module that is the direct sum of n'th tensor product,
        return TensorStarElement(dict(expanded), parent)
    else:
        TensorStarElement(dict(expanded), parent)
     
def TensorDGAlgebra(Tensor, DGAlgebra):
    '''Tensor Product of DGAlgebra(Differential Graded Algebra) is a DG Algebra'''
    def diff(self,gen):
        return E0.accumulate([
                expandTensor(gen[:i]+(gen[i].diff(),)+gen[i+1:], self)
                for i in range(len(gen))]) # CB and ask arkam multiplication? # F2?
    
    def multiply(self, gen1, gen2):
        if not isinstance(gen1, TensorGenerator) or gen1.parent!= self:
            return NotImplemented
        if not isinstance(gen2, TensorGenerator) or gen2.parent!= self:
            return NotImplemented
        return expandTensor(tuple([gen1[i]*gen2[i] for i in range(len(self))]),
                            self)
    
    def getGenerators(self):
        ''' Return the set of generators. Use product of sets of generators
        of components. Only implemented for tensors of two algebras''' #CB only 
        if len(self)!=2:
            return NotImplemented
        gens1 = self[0].getGenerators()
        gens2 = self[1].getGenerators()
        result = []
        for gen1 in gens1:
            for gen2 in gens2:
                result.append(TensorGenerator((gen1, gen2),self))
        return result
        
def TensorIdempotent(tuple): #CB and ask akram
    '''Serves as idempotent to a tensor product of algebras.'''
    def toAlgElt(self, parent):
        ''' get the algebra element corresponding to this idempotent'''
        assert isinstance(parent, TensorDGAlgebra)
        return TensorGenerator(tuple([self[i].toAlgElt(parent[i])
                                      for i in range(len(self))]), parent) 
    
def simplifyComplex(arrows, default_coeff = E0, find_homology_basis = False,
                    cancellation_constraint = None):
    '''Simplify complex using the cancellation lemma.'''
    
    '''`arrows` from simplify() function class ChainComplex. It represents the complex to be simplified. It is a dictionary object
    with generators as keys, and values are dictionary objects mapping generators to coefficients. 
    our algebra elements. For example `arrows[x][y] = coeff``` means there 
    exists an arrow from x to y with coefficient ``coeff`.
    
     The simplification is done via cancellation lemma: i.e For any arrow
     from x to y, with invertible cofficient ``coeff`` ( algebra
     element for our case, the generators x and y can be cancelled as follows:
    remove x,y and all arrows entering or leaving these two generators. 
    
    For each arrow  a- > y with coeff c1,
    each arrow x -> b with  coeff c2, in the previous complex, add an arrow
    from a to b with coefficient c_1 * coeff^-1 * c_2. The new coefficient is 
    added onto the previous one if an arrow already exists from a t o b, possibly 
    cancelling the previous arrow.
    
    `default coefficient` specifies the value to be used for zero coefficients. 
    It should be 0 if the coefficients are numbers, E0 if type Element.
    
    `find_homology_basis` specifies whether to keep track of identity of the generators.
    
    `cancellation_constraint` is a function taking two generators as
    arguments, and returns a boolean stating whether this pair of generators may
    be cancelled. If None, then any pair of generators can be cancelled. This is
    usually used to specify some condition on filtrations of generators.
    
    The technique to find optimalways to locate arrows to cancel such that maximum number
    of cancels are cancelled are as follows:
        `rev_arrows` dictionary object that for each generator y, keeps the list of 
        arrows going into y
        
        Roughly use Priority Queue - where for each arrow x to y, compute its degree 
        is defined to be (|A|-1)(|B|-1)， 
        |A| is the set of arrows going into y, 
        |B| is the set of arrows coming from x. 
        Note that the degree of an arrow may change as we apply cancellation lemma
     '''
    # Produce rev_arrows  ( basically, reverse arrows)
    rev_arrows = dict()
    for x in arrows: # print generators
        rev_arrows[x] = dict()
    for x in arrows:
        for y in arrows[x]:
            if isinstance(arrows[x][y],Element):
                rev_arrows[y][x] = arrows[x][y].copy()
            else:
                rev_arrows[y][x] = arrows[x][y]
    
    cancel_list = []
    def tryAddEdge(x,y):
        ''' If the arrow from x to y can be cancelled, then add it to cancel_list
        (a heap object), along with its degree'''
        if cancellation_constraint is not None:
            if not cancellation_constraint(x,y):
                return
        coeff = arrows[x][y]
        if coeff.invertible():
            cur_degree = (len(arrows[x]) - 1) * (len(rev_arrows[y])-1)
            heapq.heappush(cancel_list,(cur_degree, x, y)) 
            #Comment : heappush will sort by the first element of the tuple:
        
    # Make an initial list of cancellable arrows
    for x in arrows:
        for y in arrows[x]:
            tryAddEdge(x,y)
    
    # Initialize prev_meaning attribute
    if find_homology_basis:
        for x in arrows:
            x.prev_meaning = 1 * x
    
    def cancelEdge(x,y):
        ''' Cancel the edge from x to y'''
        coeff = arrows[x][y]
        assert coeff.invertible()
        inv_coeff = coeff.inverse()
        
        #List of edges going into y ( other than that from x and y)
        alist = [(term, coeff) for term, coeff in rev_arrows[y].items() if term not in (x,y)]
        #List of edges going from x ( other than that going to x and y)
        blist = [(term, coeff) for term, coeff in arrows[x].items() if term not in (x,y)]
        
        #Remove all edges going into / out of x or y
        for term in arrows[x]: # arrows out of x
            if term not in (x,y):
                del rev_arrows[term][x]
        for term in arrows[y]:# arrows out of y
            if term not in (x,y):
                del rev_arrows[term][y]
        for term in rev_arrows[x]: # arrows into x
            if term not in (x,y):
                del arrows[term][x]
        for term in rev_arrows[y]:
            if term not in (x,y):
                del arrows[term][y]
            
        # Remove x and y
        del arrows[x], arrows[y], rev_arrows[x], rev_arrows[y]
        
        # Add arrows from alist to blist
        for a, c1 in alist:
            c1_inv = c1 * inv_coeff
            for b,c2 in blist:
                new_coeff = c1_inv * c2
                if b not in arrows[a]:
                    arrows[a][b] = default_coeff
                    rev_arrows[b][a] = default_coeff
                arrows[a][b] += new_coeff
                rev_arrows[b][a] += new_coeff
                # if becomes zero
                if arrows[a][b] == default_coeff: # differnt from bzhan cb
                    del arrows[a][b], rev_arrows[b][a]
                else:
                    tryAddEdge(a,b)
       
        #update prev_meaning
        if find_homology_basis:
            for a, c1 in alist:
                a.prev_meaning += (c1 * inv_coeff) * x.prev_meaning
        
    # Main Loop: Try to cancel each edge in priority queue:
    while cancel_list:
        degree, x, y = heapq.heappop(cancel_list)
        if x in arrows and y in arrows[x] and arrows[x][y].invertible():
            new_degree = (len(arrows[x])-1)*(len(rev_arrows[y])-1)
            if new_degree > degree * 2: 
                tryAddEdge(x,y) # ask akram CB useless
            else:
                cancelEdge(x,y)
    return arrows
        
class TensorStarGenerator(Generator, tuple): # ask akram if necessary
    '''Represents a generator of the tensor star algebra - a tuple (possibly
    with zero components) of elements in the same algebra.'''
    pass

class TensorStarElement(Element):
    ''' Represents an element of the tensor star algebra.'''
    pass

class TensorStar(FreeModule): # CB ask Akram
    '''Represents a free module that is the direct sum of n'th tensor product,
    over n >= 0, of some free module A. So each generator is a (possibly empty)
    sequence of generators of A.'''
    pass
    
def findRankOverF2(num_row, num_col, entries): # CB ask Akram - necessary
    '''Find rank of a matrix over F2 with the given number of rows and columns.
    entries is a list of pairs (i, j) with 0 <= i < num_row and 0 <= j < num_col
    specifying where the matrix has 1's.'''
