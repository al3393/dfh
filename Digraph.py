'''Handles tasks related to directed graphs.'''

from algebra import Generator, SimpleChainComplex
from algebra import E0
from dstructure import DGenerator, SimpleDStructure
from ddstructure import DDGenerator, SimpleDDStructure
from dastructure import DATensorDGenerator, DATensorDDGenerator, \
    SimpleDAGenerator, SimpleDAStructure
from grading import GeneralGradingSet, GeneralGradingSetElement
from identityaa import homotopyMap
from pmc import Strands, StrandDiagram
from utility import memorize
from utility import F2

class DiGraph:
    ''' Interface for a general directed graph.'''    
    def get_Edges(self, node):
        '''Returns the list of outward edges from a node.'''
        raise NotImplementedError("Get list of outward edges not implemented.")
    
class DiGraphNode():
    '''A general node in a digraph.''' 
    def __init__(self, parent):
        ''' Endow each node a parent. '''
        self.parent = parent

class DiGraphEdge():
    ''' A general edge in a digraph.'''
    def __init__(self, source, target):
        ''' Every edge needs a source and a target (nodes from a same digraph) '''
        assert source.parent == target.parent
        self.source = source
        self.target = target

class ConcreteDiGraph(DiGraph):
    ''' Represents a directed graph as a list of nodes and list of edges.'''
    def __init__(self):
        ''' Initialize an empty concrete digraph.'''
        self.nodes = []
        self.edges = dict() # dictionary with key: nodes , values: list of edges
    
    def __str__(self):
        result = "Directed graph with %d nodes." %len(self.nodes)
        for node, outs in self.edges.item():
            node_str = "Node %s with outward edges: \n" % str(node)
            node_str += "\n".join([str(edge) for edge in outs])
            result += ("\n" + node_str)
        return result

    def addNode(self, node):
        ''' Add a node.'''
        assert node.parent == self
        self.nodes.append(node)
        assert node not in self.edges 
        self.edges[node] = []
        
    def addEdge(self, edge):
        '''Add an edge between two nodes in this directed graph. '''
        self.edges[edge.source].append(edge)
    
    def getNodes(self):
        ''' Returns the nodes of this directed graph.'''
        return self.nodes
    
    def getOutEges(self,node):
        ''' Returns the edges of the `node`.'''

class TypeDGraphNode(DiGraphNode):
    ''' A node in a type D graph. Stores the generator of type D Structure 
    corresponding to this node'''
    
    def __init__(self, parent, dgen):
        ''' Specifies the `type D structure` corresponding to this node( generator'''
        DiGraphNode.__init__(self, parent)
        self.dgen = dgen
        self.idem = self.dgen.idem # cb and remove / modify
    
    def __str__(self):
        return str(self.dgen)
    
    def __repr__(self):
        return str(self)
    
    def __eq__(self, other):
        return self.dgen == other.dgen
    
    def __ne__(self,other):
        return not (self == other)
    
    def __hash__(self):
        return hash((self.dgen, "TypeDGraphNode"))
        
class TypeDGraphEdge(DiGraphEdge):
    '''An edge in a type D graph. An edge corresponds to a type D operation. 
    Stores the algebra coefficient.'''
    
    def __init__(self, source, target, coeff):
        ''' Specifies the algebra coefficient.'''
        DiGraphEdge.__init__(self, source, target)
        self.coeff = coeff
    
    def __str__(self):
        return "Edge from %s to %s with coeff %s" % \
        (str(self.source), str(self.target), str(self.coeff))
    
    def __repr__(self):
        return str(self)
    
class TypeDGraph(ConcreteDiGraph):
    '''Corresponds to a type D structure. Nodes correspond to 
    generators and edges correspond to a type D operations.'''
    
    def __init__(self,d_structure):
        '''Creates a type D graph from type D structure.'''
        ConcreteDiGraph.__init__(self)
        
        # Maintain a dictionary from a D Structure Generator to nodes in a graph
        self.graph_node = dict()
        self.algebra = d_structure.algebra
        
        # Add nodes   
        for dgen in d_structure.getGenerators():
            cur_node = TypeDGraphNode(self, dgen)
            self.addNode(cur_node)
            self.graph_node[dgen] = cur_node
        #Add Edges
        for gen_from in d_structure.getGenerators():
            for (alg_coeff, gen_to), ring_coeff in gen_from.delta().items():
                self.addEdge(TypeDGraphEdge(self.graph.node[gen_from], \
                            self.graph_node[gen_too], alg_coeff))

class UniversalDiGraphNode(DiGraphNode, tuple):
    ''' A node in the universal directed graph of an Algebra. A sequence of Algebra
    Generators.'''
    
    def __new__(cls, parent, data):
        return tuple.__new__(cls,data)
    
    def __init__(self,parent,data):
        ''' Specifies the parent Directed Graph. Data is the list of Generators.'''
        self.parent = parent
    
class UniversalDiGraph(DiGraph):
    ''' The universal digraph of an algebra ( Strand Algebra of a tangle, for us).
    Nodes correspond to ordered sequences of algebra generators. 
    From each node, the set of outward edges is the set of algebra generators,
    appending that generator onto the sequence.'''

    def __init__(self, algebra):
        ''' Create the universal digraph for the given algebra.'''
        self.algebra = algebra
    
    def get_Edges(self, gen_from):
        result = []
        for alg_gen in self.algebra.getGenerators():
            if not alg_gen.isIdempotent():
                gen_to = UniversalDiGraphNode(self,gen_from + (alg_gen,)) #? what does '+`mean
                result.append(TypeDGraphEdge(gen_from, gen_to, alg_gen))
        return result
    def getInitialNode(self):
        '''Return the starting node, consisting of the empty sequence of generators.'''
        return UniversalDigraphNode(self,tuple())


class ATensorDGenerator(Generator, tuple):
    ''' Generator of a chain complex formed by tensoring a type A 
    structure and a type D Structure. (actually, the result of tensoring I(dL) * CT(Ti) * I(dR)'''

    def __new__(cls, parent, gen_left, gen_right):
        return tuple.__new__(cls, (gen_left, gen_right))
    
    def __init__(self, parent, gen_left, gen_right):
        ''' Specifies generators on the two sides of the tensor.''' # Both are type D generators? ( CB and reread)
        
        Generator.__init__(self, parent)


def computeATensorD(dstr1, dstr2):
    ''' Compute the tensor product dstr1 * CT(Ti) * dstr2. '''
    

def computeDATensorD(dastr1, dstr2):
    ''' Compute the tensor product ddstr1 * CT(Ti) * ddstr2.
    ddstr1 is a left '''

def computeATensorDA(dstr1, ddstr2):
    pass
                                       
        