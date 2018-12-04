import operator
from copy import deepcopy
#from algebra import DGAlgebra, Element, Generator, Tensor, TensorGenerator, 
import statistics as st
from utility import orientation,orientation_i, complement,generate_subset,\
 doescross, doescross_simple, intersections, simple_intersections, F2
from Algebra import DGAlgebra, Element, Generator, Tensor, TensorGenerator
from Algebra import E0
from utility import get_domain, get_range, generate_bijections, combinations, \
                get_domain_dict, get_range_dict, doescross_simple_rc,in_between_list, \
                reorganize_sign, reorganize_sign_2,dict_shift, get_start_dict, get_end_dict,\
                dict_shift_double, doescross_bool, mod_between, mod_helper, replace_sd_1, replace_sd_2,generate_bijections_same,\
                generate_bijections_3
'''Elementary tangles and its algebras'''

class TANGLE:
    '''represents an elementary tangle that comes after applying braid word to
    a certain tangle.'''
    def __init__(self, pairs):
        '''Creates a tangle from a dictionary object (`pairs`)of matched pairs.'''
    # `pairs` is a dictionary object which contains each orange tangle
    #in `pairs`; orientation is `key` to `value`
    
    # Note: 
    # if a the tangle is at the boundary, then we will have a dummy pair 
    # if only_right_half = =True i.e leftmost cup, then l_pairs = {}
    # if only_right_half == True i.e leftmost cup, then r_pairs = {}
               
        self.pairs = pairs
        self.n = len(pairs) # assume unique pairs 
        self.i_minus, self.i_mid, self.i_plus, self.only_left_half,self.only_right_half =  self.get_x_coord(self.pairs)  
        self.l_pairs, self.r_pairs = self.split_pairs(self.i_minus,self.i_mid, self.i_plus, self.pairs)
        self.left_boundary = self.getleft()
        self.right_boundary = self.getright()
        self.undirect_pairs_split()
        self.alpha_left = self.get_alpha_left()
        self.alpha_right = self.get_alpha_right()
        self.beta = self.get_beta()
        self.orient_left_lhalf, self.orient_right_lhalf = self.split_directions(True)
        self.orient_left_rhalf, self.orient_right_rhalf = self.split_directions(False)
                
    def __eq__(self,other):
        return self.pairs == other.pairs
    
    def __ne__(self,other):
        return not (self == other)

    def __str__(self):
        return str(self.pairs)
    
    def sd(self, data): # cb and fill
        raise NotImplementedError
    
    def split_pairs(self,a, b, c, pairs):
        '''this method takes an elementary tangle with x coordinate values,
        a, b, and c and returns two dictionaries left_pairs, right_pairs.'''   
        pairs = pairs
        l_pairs = {}
        r_pairs = {}     
        if (not(a == 0 or c == 300)):
            for key,value in pairs.items():
                if a in {key[0],value[0]}:
                    l_pairs.update({key:value}) 
                elif c in {key[0],value[0]}:
                    r_pairs.update({key:value})
        else:  # remove the dummy pair
            for key,value in pairs.items():
                if(key[0] == 0 or value[0] == 0) or (key[0] == 300 or value[0] == 300):
                    key_name = key                  
            del pairs[key_name]
            if a == 0:
                r_pairs = pairs
                l_pairs = {}       
            else:
                assert c == 300
                l_pairs = pairs
                r_pairs = {}      
        return l_pairs,r_pairs
    
    def remove_cups_caps(self):
        ''' returns two dictionary objects, l_rpairs_nc, r_rpairs_nc
        that are l_pairs and r_pairs with removed caps or caps.'''     
        temp_left = self.l_pairs
        temp_right = self.r_pairs
        l_pairs_nc = {}
        r_pairs_nc = {}            
        #left half
        l_pairs_nc = { k : temp_left[k] for k in set(temp_left) - set(self.caps) } 
        #right half
        r_pairs_nc = { k : temp_right[k] for k in set(temp_right) - set(self.cups) }
        self.l_pairs_nc = l_pairs_nc
        self.r_pairs_nc = r_pairs_nc 
        return l_pairs_nc, r_pairs_nc
    
    def add_cups_caps(self):
        ''' This method adds the cups and caps by splitting the index
        for example, if a cap is (1,5):(1,4) then it adds to 
        l_pairs_wc (1,5):(1.5,4.5) and (1.5,4.5):(1,5).'''  
      
        # remove caps or cups
        self.remove_cups_caps()        
        #add caps or cups again --------------------------------
        l_pairs_wc = self.l_pairs_nc
        r_pairs_wc= self.r_pairs_nc   
        if self.caps != {}:#left half
            if len(self.caps) != 1:
                assert TypeError("Something is wrong -- more than one caps")    
            else: # there is one cap
                for key, value in self.caps.items():            
                    # If cap
                    if orientation((key,value),True) == (-1, 1) or \
                    orientation((key,value),True) == (1, -1):                        
                        if (key[0] + 0.5) != self.i_mid:
                                raise TypeError("Something is wrong with cap length")
                        new_pt = (key[0]+ 0.5, min(key[1],value[1]) + 0.5)
                        l_pairs_wc.update({key:new_pt, new_pt: value})                        
                    else:
                        raise TypeError("Something is wrong -- this is not a cap")  
                        
        if self.cups != {}: #right half
            if len(self.cups) != 1:
                assert TypeError("Something is wrong -- more than one cups")          
            else: # there is one cup
                for key, value in self.cups.items():  
                    #if cup
                    if orientation((key,value), False) == (-1, 1) or \
                    orientation((key,value),False) == (1, -1):                  
                        if (key[0] - 0.5) != self.i_mid:
                                raise TypeError("Something is wrong with cup length")            
                        new_pt = (key[0]- 0.5, min(key[1],value[1]) + 0.5)
                        r_pairs_wc.update({key:new_pt,new_pt: value})    
                    else:
                        raise TypeError("Something is wrong -- this is not a cup")
        self.l_pairs_wc = l_pairs_wc
        self.r_pairs_wc = r_pairs_wc
        
    def split_directions(self, left_half = True):
        ''' Returns two dictionary objects, `orient_right` and 
        `orient_left` for those pairs going left and right.''' 
        self.get_cups()
        self.get_caps()
        self.remove_cups_caps() 
        self.add_cups_caps()
        
        orient_left = {}
        orient_right = {}   
        
        if left_half == True:#left_half
            pairs = self.l_pairs_wc
        else: # right_half
            pairs = self.r_pairs_wc       
            
        for key, value in pairs.items():
            if orientation((key,value), True) == (1):
                orient_right.update({key:value})
            elif orientation((key,value), True) == (-1):
                orient_left.update({key:value})
            else:
                print("Something is wrong --add_cups_and_caps() did not work properly")
        return orient_left, orient_right
                      
    def undirect_pairs_split(self): # useless cb and remove
        '''given l_pairs, r_pairs redirect the pairs such as
        it goes from left to right. fixes orientation for 
        caps or cups as well such that it always goes clockwise.
        DOES NOT SPLIT CUP OR CAP YET'''  
        self.ud_pairs_l = {}
        self.ud_pairs_r = {}   
        # left right
        for key,value in self.l_pairs.items(): 
            orient = orientation((key,value), True)     
            if orient == (1):
                self.ud_pairs_l.update({key:value})
            elif orient == (-1):
                self.ud_pairs_l.update({value:key})
            #cap
            elif orient == (1, -1):
                self.ud_pairs_l.update({value:key}) # reverse
            else:
                self.ud_pairs_l.update({key:value}) # maintain 
        # right half
        for key,value in self.r_pairs.items():
            orient = orientation((key,value), False)
            if orient == (1):
                self.ud_pairs_r.update({key:value})
            elif orient == (-1):
                self.ud_pairs_r.update({value:key})
            #cup
            elif orient == (1, -1):
                self.ud_pairs_r.update({value:key}) # reverse
            else:
                self.ud_pairs_r.update({key:value}) # maintain
        return self.ud_pairs_l, self.ud_pairs_r
    
    def undirect_pairs(self):
        split = self.undirect_pairs_split()
        return dict(split[0], **split[1])
    
    def get_x_coord(self,pairs):
        pairs = self.pairs     
        # extract i, i+0.5, i+1
        s_1 = set(val[0] for val in pairs.keys())  
        s_2 = set(val[0] for val in pairs.values())  
        s = s_1.union(s_2)
        
        assert len(s) == 3 # if it is not an elementary tangle assert error
        
        # classify whether it belongs to left half or right half
        only_left_half = False
        only_right_half = False    
        if len(s) == 3:
            i_minus = min(s)  #i
            i_mid = st.median(s) # i+0.5
            i_plus= max(s) #i + 1
        else:
            print("Something is wrong -- tangle has less than three \
                   x coordinates.")       
        if i_minus == 0: # for boundaries
            only_right_half = True #type A
        if i_plus == 300: 
            only_left_half = True # type D
            
        return i_minus, i_mid, i_plus, only_left_half, only_right_half
        
    def getleft(self):
        '''returns the `left` dictionary object containing left boundary.
            for example, for a cup dLT empty, where coordinates are {(1,0),(1,1)} 
            and dRT = {-,+}, 
            left = {}, right = {(1,0):-1, (1,1):1}.'''  
    # 'l_boundary' is an dictionary of whose key is coordinate and value is 
    #dLT sign sequence 
        left = self.l_pairs
        l_boundary = {}      
        if left == {}: # left boundary - that is only one cup
            pass
        else:
            for key, value in left.items():
                l_boundary.update(orientation_i((key,value),True))
        return l_boundary

    def getright(self):
        '''returns the `right` dictionary object containing the right boundary.'''
    # 'r_bounary' is an dictionary such that key is a coordinate and value is 
    #dRT sign sequence 
        right = self.r_pairs
        r_boundary = {}
        # right boundary - that is only one cap
        if right == {}:
            pass
        else:
            for key, value in right.items():
                r_boundary.update(orientation_i((key,value), False))
        return r_boundary
    
    def num_alpha_left(self):
        '''this method returns the size of A_i  '''
        return len(self.get_alpha_left())
    
    def num_alpha_right(self):
        ''' this method returns the size of A_{i+1}'''
        return len(self.get_alpha_right())
        
    def num_beta(self):
        ''' this method returns n the size of B_{i+1}'''
        return len(self.get_beta())
      
    def get_alpha_left(self):
        ''' returns possible alpha states:
            key : index s of alpha curve A_i
            value : coordinate'''
     # for each boundary points, make an alpha curve below( in y-axis shifted by 0.5)
        alpha = {}
        # if furthest left tangle 
        if self.left_boundary == {}:
            alpha = {0:(0,1.5)}
        else:
            assert self.left_boundary != {} #assert error if otherwise
            s_i = len(self.left_boundary) # number of tangles touching boundary
            coord = list(self.left_boundary.keys())
            coord.sort() #sort coordinates, in case it is not sorted 
            for i in range(s_i):
                x, y = coord[i]
                alpha.update({i:(x,(y-0.5))})       
            # add the top alpha curve
            x, y = coord[s_i -1]
            alpha.update({s_i: (x, y+ 0.5)})
        return alpha
    
    def get_alpha_right(self):
        ''' returns possible alpha states:
            key : index is of alpha curve A_{i+1}
            value : coordinate '''
        alpha = {}
        if self.right_boundary == {}:  # if furthest left tangle 
            alpha = {0:(300,1.5)}
        else:
            assert self.right_boundary != {} #assert error if otherwise
            s_i = len(self.right_boundary) # number of tangles touching boundary
            coord = list(self.right_boundary.keys())
            coord.sort() #sort coordinates, in case it is not sorted
            for i in range(s_i):
                x, y = coord[i]
                alpha.update({i:(x,(y-0.5))})
            x, y = coord[s_i -1] # add the top alpha curve
            alpha.update({s_i: (x, y+ 0.5)}) 
        return alpha
    
    def occupied_B(self):
        '''uses information on pairs and returns the set object of 
        y-coordinates occupied in B by the tangle, including a cup or a cap.'''
        self.occupied = set() # set of occupied y axis
        for value in list(self.ud_pairs_l.values()): #contribution from left half
            if value[0] == self.i_mid: # not cap
                self.occupied.add(value[1])
            else:  # cap
                #raise value[0] != self.i_mid
                self.occupied.add(value[1] + 0.5)
        for key in list(self.ud_pairs_r.keys()): #contribution from right half
            if key[0] == self.i_mid: # not cup 
                self.occupied.add(key[1])  
            else: # cup
                #raise key[0] != self.i_mid
                self.occupied.add(key[1] + 0.5)
        return self.occupied
    
    def get_beta(self):
        ''' returns a dictionary object of possible beta states:
            key : index s of alpha curve B_{i+1}
            value : coordinate''' 
        self.occupied_B()
        b_occupied = list(self.occupied)
        b_occupied.sort()   
        beta = {}
        s_i = len(b_occupied)   
        for i in range(s_i):
            y = b_occupied[i]
            beta.update({i:(self.i_mid, y - 0.5)})
        # add top beta curve
        y = b_occupied[s_i - 1]
        beta.update({s_i: (self.i_mid, y + 0.5)})  
        return beta
    
    def get_caps(self):
        '''fetches any caps that may exist.'''
        left = self.l_pairs
        self.caps = {}
        for key,value in left.items():
            if key[0] == value[0]:
                self.caps.update({key:value})          
        return self.caps
    
    def get_cups(self):
        '''fetches any cups that may exist.'''
        right = self.r_pairs
        self.cups = {}
        for key,value in right.items():
            if key[0] == value[0]:
                self.cups.update({key:value})
        return self.cups
    
    def complement_idem(self,data, is_left):
        '''Given an set of occupied state alpha arcs (left or right), it 
        returns the 'unoccupied states' as tuple object `data` '''     
        if is_left == True:
            alpha = tuple(sorted(i for i in range(self.num_alpha_left())))
            return complement(alpha,data)
        else:
            alpha = tuple(sorted(i for i in range(self.num_alpha_right())))
            return complement(alpha,data)
             
    def idem(self,data):# CB and remove
        raise NotImplementedError
#    def grading(self, maslov, alexander): #cb 
#        grading_group = GradingGroup(self)
#        return GradingElement(grading_group, maslov, alexander)
    
    def big_gr(self,maslov, alexander):# cb and remove 
        raise NotImplementedError
        
    def is_simple_tangle(self):
        '''adds true to the variable `is_simple`, when the tangle is a simple tangle' used for
        Simple_Strand'''
        self.is_simple = True
    
    def getIdempotents(self, is_left): 
        '''Get the list of idempotents depending on variable `is_left`.
        If `is_left` is true, takes into account of left boundary of this tangle''' 
        idems = [] 
        if is_left == True:
            n = self.num_alpha_left() - 1 
        else:
            n = self.num_alpha_right() - 1 
        lst = generate_bijections(n)
        for perm in lst:
            idems.append(Idempotent(self, perm, is_left))      
        return idems
    
    def getAlgebra(self, is_left, parent): #CB and check
        '''Get the list of generators of the strand algebra, depending on variable `is_left`.
        If `is_left` is true, takes into account of left boundary of this tangle
        `parent` is the Strand Algebra of corresponding `is_left` side '''

        if self.only_left_half and is_left == False:
            return TypeError("TypeError: This is a type D Tangle .")
        if self.only_right_half and is_left == True:
            return TypeError("TypeError: This is a type A Tangle.")
        alg_gen = []
        if is_left:
            n = len(self.alpha_left)-1
        else:
            n = len(self.alpha_right) -1
        bij = generate_bijections_same(n)
        for bijection in bij:
            alg_gen.append(Simple_Strand(parent, is_left, tuple(bijection)))
        return alg_gen

    def getStrandDiagrams(self,p_module): # cb and check 
        ''' Get the list of generators of CT(Ti) ( the parent module). '''
        gen = []
        n = len(self.alpha_left)-1
        m = len(self.beta) - 1
        k = len(self.alpha_right) -1
        bij = generate_bijections_3(n,m,k)
        for bijection in bij:
#            print("bijection is: {0}".format(bijection))
            gen.append(StrandDiagram(p_module, Strands(self,bijection)))
        return gen
    
class Idempotent(tuple):
    '''Represents an idempotent in a certain tangle. Stored as a tuple of occupied pairs.'''
    def __new__(cls, p_tangle, data, is_left = True):
        '''If strands are going from 1 - > 1, 3 -> 3, then data will consist of (1,3)''' 
        return tuple.__new__(cls, tuple(sorted(data)))
    
    def __init__(self, p_tangle, data, is_left):
        self.p_tangle = p_tangle
        self.is_left = is_left
        
        if is_left and self.p_tangle.only_right_half:
            raise TypeError("Left empty boundary")
        if not is_left and self.p_tangle.only_left_half:
            raise TypeError("Right empty boundary")
            
        if self.is_left == True:
            self.sign_seq = self.p_tangle.left_boundary
        else:
            self.sign_seq = self.p_tangle.right_boundary
        self.tangle = TANGLE(self.get_tangle_pairs())
        self.tangle.is_simple_tangle()
        self.pairs = self.get_tup()
    
    def __eq__ (self, other): # cb and change - idempotents tangle when moving
        if isinstance(other, Idempotent):
            return self.sign_seq == other.sign_seq and tuple.__eq__(self,other)
        elif isinstance(other, tuple):
            return tuple.__eq__(self,other)
        else:
            return False
        
    def __ne__(self,other):
        return not (self == other)
    
    def __hash__(self):
        return hash((self.tangle, tuple(self), "Idempotent"))
    
    def __str__(self):
        return repr(self)
    
    def __repr__(self):
        
        if self.is_left == True:
            return "(%s)" % " , ".join(str(i) + ":" + str(self.p_tangle.alpha_left[i]) for i in self)
        else:       
            return "(%s)" % " , ".join(str(i) + ":" + str(self.p_tangle.alpha_right[i]) for i in self)
    
    def get_tangle_pairs(self):
        ''' returns the dictionary object that is needed to make a new tangle.'''   
        pairs = {}
        for coord,sign in self.sign_seq.items():
            if self.is_left:
                np1 = (coord[0] - 0.5, coord[1])
                np2 = (coord[0] - 1, coord[1])
                if sign == -1:
                    pairs.update({coord: np1, np1: np2})
                else:
                    if sign != 1:
                        raise ValueError("Something is wrong with the sign")
                    pairs.update({np2:np1, np1:coord})
            else:
                np1 = (coord[0] + 0.5, coord[1])
                np2 = (coord[0] + 1, coord[1])
                if sign == -1:
                    pairs.update({np2:np1, np1:coord})      
                else:
                    if sign != 1:
                        raise ValueError("Something is wrong with the sign")
                    pairs.update({coord: np1, np1: np2})
        return pairs
   
    def get_tup(self):
        return tuple((k,k) for k in self)
            
    def comp_idem(self):
        ''' Get the complementary idempotent in the same tangle on the same side.'''
        return Idempotent(self.p_tangle, self.p_tangle.complement_idem(self, self.is_left), \
                                                self.is_left) 
    def toAlgElt(self,parent):
        ''' Get the strand algebra element corresponding to this idempotent
        in the specified `parent` strand algebra'''
        if parent.tangle != self.p_tangle:
            raise TypeError("Something is wrong -- tangle T_i don't match")
        else:
            return Simple_Strand(parent,self.is_left ,self.get_tup())
              
class Strands(tuple):
    ''' Represents a (fixed) list of strands in a certain tangle. Stored as
    a tuple of tuple of pairs.
    For ex) data will be of form: (((2,1),(4,3),(5,4)),((2,1),(5,4)))
    where first tuple contains the strands on the left, first index alpha 
    where second tuple contains the strands on the right,second index alpha'''  
    def __new__(cls, tangle, data):
        data_0 = tuple(sorted(data[0]))
        data_1 = tuple(sorted(data[1]))
        return tuple.__new__(cls, tuple([data_0, data_1]))

    def __init__(self, tangle, data):
#        data_0 = tuple(sorted(data[0]))
#        data_1 = tuple(sorted(data[1]))
        self.tangle = tangle
        self.data = data
        self.convert() # cb and remove later
        self.converted  = self.left_converted.copy() #dictionary format
        self.converted.update(self.right_converted)
        self.left_crossings = self.strandCrossing(True) #list of crossings on the left
        self.right_crossings = self.strandCrossing(False)        
     
    def __eq__(self,other):
        if isinstance(other, Strands):
            return self.tangle == other.tangle and tuple.__eq__(self,other)
        else:
            return False
    
    def __ne__(self,other):
        return not (self == other)
    
    def __hash__(self):
        return hash((self.tangle, tuple(self), "Strands"))
    
    def __str__(self):
        string = "\n Left: "
        string += "(%s)" % ",".join("%s->%s" % (s, t) for s, t in self[0]) 
        string += "Right: "    
        string += "(%s)" % ",".join("%s->%s" % (s, t) for s, t in self[1]) 
        return string
    
    def __repr__(self):
        return str(self)
  
    def occupied_left_alpha(self):
        '''returns a tuple of occupied left alphas curves'''
        occupied = []
        for strand in self[0]:
            occupied.append(strand[0])
        return tuple(sorted(occupied))
    
    def occupied_right_alpha(self):
        '''returns a tuple of occupied left alphas curves'''    
        occupied = []
        for strand in self[1]:
            occupied.append(strand[1])
        return tuple(sorted(occupied))
    
    def beta(self, is_left = True):
        '''returns a tuple of occupied beta curves, 
        if `is_left` then returns a tuple of occupied beta
        on the left. '''    
        occupied = []
        if is_left:
            for strand in self[0]:
                occupied.append(strand[1])
        else:
            for strand in self[1]:
                occupied.append(strand[0])
        return tuple(sorted(occupied))
    
    def leftCompatible(self,sa_gen):
        '''Tests whether this set of strands is left compatible with a given
        A(-dLT) algebra element.'''
    
        if (not isinstance(sa_gen, Simple_Strand)) and sa_gen.is_left == True:
            raise TypeError("Algebra Element and Strand not Compatible.")  
        return sa_gen.t == get_domain(self.get_left_idem().pairs())
        
    def rightCompatible(self,sa_gen):
         '''Tests whether this set of strands is compatible with a given 
         A(dRT) algebra element.'''
         if (not isinstance(sa_gen, Simple_Strand)) and sa_gen.is_left == False:
            raise TypeError("Algebra Element and Strand not Compatible.")
         return sa_gen.s == get_range(self.get_right_idem().pairs())

    def leftCompatible_idem(self,idem):
        ''' Test whether this set of strands is compatible with a given left
        idempotent'''
        raise isinstance(idem, Idempotent) and idem.is_left == True
        return self.occupied_left_alpha() == idem.comp_idem() 
    
    def rightCompatible_idem(self,idem):
        ''' Test whether this set of strands is compatible with a given right
        idempotent'''  
        raise isinstance(idem, Idempotent) and idem.is_left == False
        return self.occupied_right_alpha()== idem

    def idemCompatible(self, left_idem, right_idem):
        ''' Tests whether this set of strands is compatible with the given 
        idempotents on the two sides. '''
        
        return self.leftCompatible_idem(left_idem) and \
        self.rightCompatible_idem(right_idem) #cb
    
    def get_left_idem(self): # For tangles # needs modification
        '''Find the left_idem given the strand information.'''
        occupied_l = self.occupied_left_alpha()
        return Idempotent(self.tangle, self.tangle.complement_idem(occupied_l,True), True)
      
    def get_right_idem(self):# For strand diagram # needs modification
        '''Find the left_idem given the strand information.'''
        occupied_r = self.occupied_right_alpha()
        return Idempotent(self.tangle, occupied_r, False)
    
    def crossings(self,option, left_half = False, right_half = False):
        ''' Gives the number of crossings between strands and/or tangles
        using the given options used for
        Alexander and Maslov gradings.
        For options refer to crossingtypes.jpg '''
        #adapted from Strand Diagram Class
        num = 0
        if option == 1:
            # tangle left & strand
            if left_half == True: ## left half:
                dict1 = self.tangle.orient_left_lhalf
                dict2 = self.left_converted
                num+= simple_intersections(dict1,dict2,False)
                
            if right_half == True: ## right half:
                dict1 = self.tangle.orient_left_rhalf
                dict2 = self.right_converted
                num+= simple_intersections(dict1,dict2,False)       
        if option == 2:
            #tangle right & strand
            if left_half == True: ## left half:
                dict1 = self.tangle.orient_right_lhalf
                dict2 = self.left_converted
                num+= simple_intersections(dict1,dict2,False)  
                
            if right_half == True: ## right half:
                dict1 = self.tangle.orient_right_rhalf
                dict2 = self.right_converted
                num+= simple_intersections(dict1,dict2,False)     
        return num
    
    def strandCrossing(self, left_half = True): 
        ''' Returns the list of crossings between moving strands
        either left_half or right_half depending on `left_half`
        For example, if in the left: crossing pairs are (4,2) and (1,3)
        and if in the right: crossing pairs are (4,1) and (1,5)
        if left_half = True: 
            returns [((4,2),(1,3)),...]
        
        if left_half = False:
            returns [((4,1), (1,5)),...]
        '''
        if left_half:  #left half
            left_half = []
            left_strands = self.data[0] # tuple of left strands         
            l_num = len(left_strands) # of left strands
            print("left_strands")
            print(left_strands)
            
            combinations = generate_subset(l_num - 1 , 2)
            for combo in combinations:
                if doescross_simple(left_strands[combo[0]],left_strands[combo[1]]) == True:
                    left_half.append((left_strands[combo[0]],left_strands[combo[1]]))
            return left_half
        else:  #right half
            right_half = []
            right_strands = self.data[1] # tuple of right strands
            print("right_strands")
            print(right_strands)
            
            r_num = len(right_strands) # of right strands
            combinations = generate_subset(r_num - 1 , 2)
            for combo in combinations:
                if doescross_simple(right_strands[combo[0]],right_strands[combo[1]]) == True:
                    right_half.append((right_strands[combo[0]],right_strands[combo[1]]))
            return right_half
  
    def numCrossing(self, left_half = True):
        if left_half:
            return len(self.left_crossings)
        else:
            return len(self.right_crossings)
    
    def strandCrossing_coord(self,left_half = True):
        '''like strandCrossing, but returns coordinates instead of alpha/beta
        states. For example, {((2,3.5),(2.5,4.5)):((2,5.5),(2.5,1))}'''
        if left_half:
            raw = self.left_crossings
        else:
            raw = self.right_crossings
        converted = {}      
        for k, v in raw:
            new_k = self.convert_s_to_coord(k, left_half)
            new_v = self.convert_s_to_coord(v, left_half)
            converted.update({new_k: new_v})
        return converted
    
    def convert_s_to_coord(self, tup, is_left = True): # check remove
        ''' Convert Strand to Coordinates. That is,
        given a tuple pair say (2,3) it gives tuple of corresponding coordinates
        in tuples. if is_left == True,
        if 2 refers to (2,3.5) and 
        if 3 refers to (2.5,4.5) in actual coordinates,
        then returns (((2,3.5),(2.5,4.5))'''    
        if is_left:
            left = self.tangle.alpha_left
            right = self.tangle.beta
        else:
            left = self.tangle.beta
            right = self.tangle.alpha_right       
        return (left[tup[0]],right[tup[1]])

    def convert_di(self):
        ''' Like convert() below, but double indexed:
        For example if (1,2) in left half, it would return {(1,2):((1,0.5),(1.5,1.5))}
        with natural orientation from left to right'''   
        left_converted = {}
        right_converted ={}
        for strand in self[0]:#left half
            start = self.tangle.alpha_left[strand[0]]
            end = self.tangle.beta[strand[1]]
            left_converted.update({strand :(start,end)})
        for strand in self[1]: #right half
            start = self.tangle.beta[strand[0]]
            end = self.tangle.alpha_right[strand[1]]
            right_converted.update({strand:(start,end)})     
        self.left_converted_di = left_converted
        self.right_converted_di = right_converted
        
    def get_strand_index(self, coord,is_left):
        ''' receives the coord of a strand, and returns it in a alpha/beta index
        format. If coord = ((1,0.5),(1.5,1.5)), returns say (1,3)'''
        self.convert_di()
        if is_left:
            for k,v in self.left_converted_di.items():
                if v == coord:
                    return k
        else:
            for k,v in self.right_converted_di.items():
                if v == coord:
                    return k
            
    def convert(self):
        ''' Converts a given strand in a tuple format, to a dictionary format
        using the coordinates stored in its parent tangle.
        For example if (1,2) in left half, it would return {(1,0.5):(1.5,1.5)}
        with natural orientation from left to right'''     
        left_converted = {}
        right_converted ={}
        for strand in self.data[0]: #left half
            start = self.tangle.alpha_left[strand[0]]
            end = self.tangle.beta[strand[1]]
            left_converted.update({start:end})
        for strand in self.data[1]:#right half
            start = self.tangle.beta[strand[0]]
            end = self.tangle.alpha_right[strand[1]]
            right_converted.update({start:end})
        self.left_converted = left_converted
        self.right_converted = right_converted
    
    def convert_dict(self,is_left = True):
        ''' Converts a given strand in a tuple format to a dictionary format
        using the coordinates stored in its parent tangle.
        For example if is_left is true, and left half of the strand constists of 
        (1,2),(3,4) in left half, it become {1:2, 3:4}'''
        left_converted = {}
        right_converted ={}
        if is_left: #left half
            for strand in self[0]:
                left_converted.update({strand[0]:strand[1]})
            return left_converted
        else:  #right half
            for strand in self[1]:
                right_converted.update({strand[0]:strand[1]}) 
            return right_converted

class StrandDiagram(Generator):
    ''' Represents a strand diagram of a tangle T_i
       i.e., a generator of CT(T_i).''' 
    # Note: need to check len(self) == len(self.tangle.beta)
    def __init__(self, parent, strands, left_idem = None, right_idem = None):
        '''Specifies Tangle, parent, and strands as a list of pairs'''
        # what does parent here mean cb and ask Akram
        
        Generator.__init__(self,parent)
        self.tangle = parent.tangle
        self.strands = strands
        if not isinstance(self.strands, Strands):
            self.strands = Strands(self.tangle, self.strands)
        if (len(self.strands[0]) + len(self.strands[1])) != len(self.tangle.beta):
            raise TypeError("The beta curves are not filled.")
        if not self.tangle.only_right_half: # Calculate left idempotent
            self.left_idem = self.strands.get_left_idem()  
        if not self.tangle.only_left_half: # Calculate right idempotent
            self.right_idem = self.strands.get_right_idem()      
#        print("+++++++++++")
#        print("\n")
#        print(self.strands)
#        print("TANGLE INFO %%%%%%%%%%")
        self.maslov = self.maslov()
        self.alexander = self.alexander()
    
    def __eq__(self, other):
        return self.parent == other.parent and self.strands == other.strands
    
    def __ne__ (self, other):
        return not (self == other)
        
    def __hash__(self):
        return hash((self.parent, self.strands))
    
    def __str__(self):
        return str(self.strands)
    
    def __repr__(self):
        return str(self)
   
    def numCrossing(self, left_half = True): # cb and remove - already in strands class
        ''' Returns the number of crossings between moving strands
        either left_half or right_half depending on `left_half'''
        if left_half ==True: #left half
            left_half = 0
            left_strands = self.strands[0] # tuple of left srands
            l_num = len(left_strands) # of left strands
            combin = generate_subset(l_num - 1 , 2)
            for combo in combin:
                if doescross_simple(left_strands[combo[0]],left_strands[combo[1]]) == True:
                    left_half += 1
            return left_half
        else: #right half  
            right_half = 0
            right_strands = self.strands[1] # tuple of right srands          
            r_num = len(right_strands) # of right strands
            combin = generate_subset(r_num - 1 , 2)
            for combo in combin:
                if doescross_simple(right_strands[combo[0]],right_strands[combo[1]]) == True:
                    right_half += 1
            return right_half
    
    def crossings(self,option, left_half = False, right_half = False):
        ''' Gives the number of crossings between strands and/or tangles
        using the given options used for Alexander and Maslov gradings.
        For options refer to crossingtypes.jpg '''     
        num = 0
        if option == 1:
            # tangle left & strand
            if left_half == True: ## left half:
                dict1 = self.tangle.orient_left_lhalf
                dict2 = self.strands.left_converted
                num+= simple_intersections(dict1,dict2,False)
            
            if right_half == True: ## right half:
                dict1 = self.tangle.orient_left_rhalf
                dict2 = self.strands.right_converted
                num+= simple_intersections(dict1,dict2,False)           
        if option == 2:
            #tangle right & strand
            if left_half == True: ## left half:
                dict1 = self.tangle.orient_right_lhalf
                dict2 = self.strands.left_converted
                num+= simple_intersections(dict1,dict2,False)
            
            if right_half == True: ## right half:
                dict1 = self.tangle.orient_right_rhalf
                dict2 = self.strands.right_converted
                num+= simple_intersections(dict1,dict2,False) 
                
        if option == 3:
            # tangle right & tangle right       
            if left_half == True: ## left half:
                dict1 = self.tangle.orient_right_lhalf
                num += simple_intersections(dict1,dict1,True)
            
            if right_half == True: ## right half:
                dict2 = self.tangle.orient_right_rhalf
                num += simple_intersections(dict2,dict2,True)      
        if option == 4:
            # tangle left & tangle_left
            if left_half == True:# left half:
                dict1 = self.tangle.orient_left_lhalf
                num += simple_intersections(dict1,dict1,True)
           
            if right_half == True:# right half:
               dict2 = self.tangle.orient_left_rhalf
               num += simple_intersections(dict2,dict2,True)     
        if option == 5:
            # tangle oriented left on entire tangle for last term in alexander
            num += len(self.tangle.orient_left_rhalf) + \
                    len(self.tangle.orient_left_lhalf)            
#            for value in self.tangle.orient_left_rhalf.values():
#                if value in self.tangle.orient_left_lhalf.keys():
#                    num += 1      
        if option == 6:
            # tangle_left on left half tangle for maslov (xi-)
            num = len(self.tangle.orient_left_lhalf) 
        return num
                    
    def maslov(self):
        '''Maslov grading for a generator of Strand Diagram'''
        maslov = 0 
#        DEBUGGING CODE
#        print("left_half (x_i -)")
#        print("Type 0 crossing:{0}".format(self.strands.numCrossing(True)))
#        print("Type 1 crossing:{0}".format(self.crossings(1, True, False)))
#        print("Type 4 crossing:{0}".format(self.crossings(4, True, False)))
#        print("Type 6 crossing:{0}".format(self.crossings(6, True, False)))
#        
#        print("right_half((x_i +)")
#        print("Type 0 crossing:{0}".format(self.strands.numCrossing(False)))
#        print("Type 2 crossing:{0}".format(self.crossings(2, False, True)))
#        print("Type 3 crossing:{0}".format(self.crossings(3, False, True)))
#        
        #left_half (x_i -)
        maslov += - self.strands.numCrossing(True) + self.crossings(1, True, False) - \
                    self.crossings(4, True, False) - self.crossings(6, True, False)  
        
        #right_half((x_i +)
        maslov += self.strands.numCrossing(False) - self.crossings(2, False, True) + \
                    self.crossings(3, False, True)
        return maslov
    
    def alexander(self):
        '''Alexander grading for a generator of Strand Diagram'''
        alex = 0
        
        #DEBUGGING CODE
#        print("ALEXANDER:")
#        print("Type 1 crossing:{0}".format(self.crossings(1, True, True)))
#        print("Type 2 crossing:{0}".format(self.crossings(2, True, True)))
#        print("Type 3 crossing:{0}".format(self.crossings(3, True, True)))
#        print("Type 4 crossing:{0}".format(self.crossings(4, True, True)))
#        print("Type 4 crossing:{0}".format(self.crossings(5, True, True)))
        
        alex += self.crossings(1, True, True) - self.crossings(2, True, True) + \
                        self.crossings(3, True, True) - self.crossings(4, True, True) \
                                    - self.crossings(5,True, True)
        return alex * 0.5 
         
    def getGrading(self):
        'return tuple (M,A) of this strand diagram'''
        return (self.maslov,self.alexander)
    
    def getLeftIdem(self):
        return self.strands.get_left_idem()
    
    def getRightIdem(self):
        return self.strands.get_right_idem()
    
    def replace(self, old, new, is_left):
        '''creates a new stranddiagram, with a new strand, with deleting pair
        `old` on the `is_left`side and replacing it with a new one.
        `old` and `new` can be an array too ''' 
        raw = deepcopy(self.strands.data)
        mod = replace_sd_1(raw, old ,new, is_left)
        return StrandDiagram(self.parent, mod)
    
    def replace_2(self, old, new):
        ''' Upgraded version for `self_2` with input as list objects 
        Old, and New is a `list` object,then do it accordingly.
        By default, the number of old and new must equal, and be in order
        they are to be replaced, respectively.'''
        
        raw = deepcopy(self.strands.data)
        mod = replace_sd_2(raw, old, new)
        return StrandDiagram(self.parent, mod)

class Simple_Strand(Generator):
    ''' Represents a generator of Strand Algebra. i.e 
    the bijections from S C [n] to T C [n]'''

    def __init__(self, parent, is_left, pairs):
        '''Specifies tangle that it takes it sign sequence from 
         parent Strand Algebra, and `pairs` as a list of pairs. Creates
         the tangle that it will use previous methods
         `pairs` is a tuple, that specifies the pairs
         if simple_strand connects 1 to 3, 3 to 4, it will contain
         ((1,3), (3,4))'''
        # parent here refers to the Strand Algebra
        Generator.__init__(self,parent)
        self.p_tangle = parent.tangle # tangle that it takes its sign from 
        self.is_left = is_left
        if is_left and self.p_tangle.only_right_half:
            raise TypeError("Left empty boundary")
        if not is_left and self.p_tangle.only_left_half:
            raise TypeError("Right empty boundary")
        if self.is_left:
            self.sign_seq = self.p_tangle.left_boundary
        else:
            self.sign_seq = self.p_tangle.right_boundary
        #check S and T in petkova paper:
        self.check_s_and_t(pairs)
        # Instantiate Base Tangle and strand: 
        self.pairs = pairs
        self.tangle = TANGLE(self.get_tangle_pairs())
        self.tangle.is_simple_tangle()
        self.strands = Strands(self.tangle, self.get_tuple_strand(self.pairs))
        self.maslov = self.maslov()
        self.alexander = self.alexander()
    
    def __eq__(self, other): # changed
        return self.sign_seq == other.sign_seq and self.pairs == other.pairs
    
    def __ne__ (self, other):
        return not (self == other)
        
    def __hash__(self):
        return hash((self.parent, self.strands))
    
    def __str__(self):
        return str(self.strands)
    
    def __repr__(self):
        return str(self)
    
    def copy(self):
        return Simple_Strand(self.parent, self.is_left, self.pairs)
    
    def replace(self,delete,add):
        '''returns a new Simple Strand generator with `delete` strand pair removed,
        and adding `add`.
            `Delete` and `add` can be a list object or tuple object.'''   
        new_pairs = []
        old_pairs = list(self.pairs)
        for pair in old_pairs:
            if pair not in delete:
                new_pairs.append(pair)      
        new_pairs = list(new_pairs) + list(add)
        new_pairs = tuple(new_pairs)
        return Simple_Strand(self.parent, self.is_left, new_pairs)
   
    def getLeftIdem(self):  # ask hsin pei CB
        return Idempotent(self.p_tangle, tuple(self.s), self.is_left)     

    def getRightIdem(self):
        return Idempotent(self.p_tangle, tuple(self.t), self.is_left)   
    
    def check_s_and_t(self,pairs):
        '''Helper method for __init__ to check whether a Strand Algebra Generator
        is well defined.'''
        self.s = get_domain(pairs) # S in petkova paper
        self.t = get_range(pairs)# T in petkova paper
        
        # making sure S, T C [n]
        if self.is_left:
            if self.p_tangle.num_alpha_left() <= max(self.s) or \
            self.p_tangle.num_alpha_left() <= max(self.t):
                raise ValueError("The strands are out of bound.")
        else:
            if self.p_tangle.num_alpha_right() <= max(self.s) or \
            self.p_tangle.num_alpha_right() <= max(self.t):
                raise ValueError("The strands are out of bound.")
        
        #Check if S and T have the same size ( i.e bijection )
        if len(self.s) != len(self.t):
            raise ValueError("Something is wrong -- S and T do not have the same size.")
            
    def get_tangle_pairs(self): # in Simple Strand
        ''' returns the dictionary object that is needed to make a new tangle.'''
        pairs = {}
        for coord,sign in self.sign_seq.items():
            if self.is_left:
                np1 = (coord[0] - 0.5, coord[1])
                np2 = (coord[0] - 1, coord[1])
                if sign == -1:
                    pairs.update({coord: np1, np1: np2})
                else:
                    if sign != 1:
                        raise ValueError("Something is wrong with the sign")
                    pairs.update({np2:np1, np1:coord})
            else:
                np1 = (coord[0] + 0.5, coord[1])
                np2 = (coord[0] + 1, coord[1])
                if sign == -1:
                    pairs.update({np2:np1, np1:coord})          
                else:
                    if sign != 1:
                        raise ValueError("Something is wrong with the sign")
                    pairs.update({coord: np1, np1: np2})             
        return pairs
    
    def get_tuple_strand(self,pairs): 
        ''' returns the tuple that is needed for strand object '''
        #if self.is_left true, then make bijections look  --- then up/down
        if self.is_left == True:
            right_half = tuple((k,k) for k in self.s)
            return tuple(( right_half, pairs))  
        # if self.is_left False, then make bijections look like up/down then ----
        else:
            left_half = tuple((k,k) for k in self.t)
        print("*********")
        print(tuple((pairs, left_half)))
        print("*********")
        return tuple((pairs, left_half))
        
    def crossings(self,option, left_half = False, right_half = False):
        ''' Gives the number of crossings between strands and/or tangles
        using the given options used for
        Alexander and Maslov gradings.
        For options refer to crossingtypes.jpg '''
        #adapted from Strand Diagram Class
        num = 0
        if option == 1:
            # tangle left & strand
            if left_half == True: ## left half:
                dict1 = self.tangle.orient_left_lhalf
                dict2 = self.strands.left_converted
                num+= simple_intersections(dict1,dict2,False)
            
            if right_half == True: ## right half:
                dict1 = self.tangle.orient_left_rhalf
                dict2 = self.strands.right_converted
                num+= simple_intersections(dict1,dict2,False)
        if option == 2:
            #tangle right & strand
            if left_half == True: ## left half:
                dict1 = self.tangle.orient_right_lhalf
                dict2 = self.strands.left_converted
                num+= simple_intersections(dict1,dict2,False)
            
            if right_half == True: ## right half:
                dict1 = self.tangle.orient_right_rhalf
                dict2 = self.strands.right_converted
                num+= simple_intersections(dict1,dict2,False)            
        return num
    
    def maslov(self): 
        '''maslov grading for a strand algebra generator'''
        return self.strands.numCrossing(True) + self.strands.numCrossing(False) - self.strands.crossings(1, True, True) 

    def alexander(self):
        '''maslov grading for a strand algebra generator'''
        return 0.5 * (self.strands.crossings(1,True,True) - self.strands.crossings(2, True, True))
    
    def isIdempotent(self):
        return self.s == self.t # self.s and self.t is tuple element, so ordering matters

    def numCrossing(self): # cb and remove - already in strands class
        ''' Returns the number of crossings between moving strands
        either left_half or right_half depending on `left_half'''
        return self.strands.numCrossing(True) + self.strands.numCrossing(False)
    
class TangleModule(DGAlgebra):
    '''CT(Ti)''' 
    def __init__(self, ring, tangle):
        '''Specifies the tangle T adn the ring'''
        DGAlgebra.__init__(self,ring)
        self.tangle = tangle
    
    def __str__(self):
        return "CT(Ti) over %s" %(str(self.tangle))
        
class StrandAlgebra(DGAlgebra): #The `parent` of  Strand Algebra
    '''Represents the strand algebra of a tangle T_i generated by 
     left or right sign sequence.'''
     
    def __init__(self, ring, tangle, is_left = True):
        '''Specifies the tangle T, and whether it takes the left or right sequence
        i.e if is_left == True, it is A(-dLT)
        if isleft == False, it is A(dRT) '''
        DGAlgebra.__init__(self,ring)
        self.tangle = tangle
        self.is_left = is_left
    
    def __str__(self):
        if self.is_left == True:
            return "Strand Algebra over %s with taking sign sequence of \
                left boundary." %(str(self.tangle))
        else:
            return "Strand Algebra over %s with taking sign sequence of \
                right boundary." %(str(self.tangle))
    
    def __eq__(self,other):
        if not isinstance(other, StrandAlgebra):
            return False
        return self.tangle == other.tangle and self.is_left == other.is_left
    
    def __ne__(self,other):
        return not (self == other)
    
    def __hash__(self):
        return hash((self.tangle, self.is_left))
     
    def getStrandDiagram(self,pairs): # CB and remove
        ''' Pairs consist of tuples ((1,3), (2,4)) etc.. i.e., Bijections of
        S to T C [n]'''
        return Simple_Strand(self, self.is_left, pairs)
    
    def diff(self, gen):  
        result = E0
        if self.ring is F2:
            for (s1, s2), dgen_term in self.diffRaw(gen):
                result += dgen_term.elt()
        else:
            return NotImplemented
        return result
    
    def diffRaw(self,gen):
        '''Returns a list of elements of the form ((s1, s2), diff_term), where
        s1 < s2 are starting points of strands in gen that crosses, and
        diff_term is a generator in gen.diff() obtained by uncrossing these two
        strands. Together they specify all terms in gen.diff(). ''' 
        lst = []
        valid_crossings = self.diff_helper_1(gen) # pairs of crossing strands that are not in figure 6 relations
        for cross in valid_crossings:
            s1 = cross[0][0]
            s2 = cross[1][0]
            a1,a2,b1,b2 = reorganize_sign(cross) # split into 4 points 
            diff_term = gen.replace(cross,((a1,b2),(b1,a2))) # replace cross with ((a1,b2),(b1,a2))
            lst.append((s1,s2, diff_term)) 
        return lst
            
    def diff_helper_1(self,gen):
        '''helper method for diffRaw:
         returns pairs of crossing strands that are not in figure 6 relations
        for cross in valid_crossings:'''   
        arr = []     
        is_left = not(gen.is_left) 
        strands_crossing = gen.strands.strandCrossing(is_left)
        gen.strands.convert_di()
        if is_left: 
            dic = gen.strands.left_converted_di
            tangles = list(gen.tangle.l_pairs.items())
            t1 = gen.tangle.i_minus
            t2 = gen.tangle.i_mid
        else:
            dic = gen.strands.right_converted_di 
            tangles = list(gen.tangle.r_pairs.items())
            t1 = gen.tangle.i_mid
            t2 = gen.tangle.i_plus
                      
        if abs(t1 - t2) >= 1:
            raise TypeError("Something is wrong- tangle width is longer than 0.5")
        
        for s_pair in strands_crossing: # for each (s1,s2)
            #print("STRAND PAIR: {0}".format(s_pair))
            a_1,a_2,b_1,b_2 = reorganize_sign_2(s_pair, dic) # cb and change  assume a1, a2 are coordinates
            inter = doescross((a_1,a_2),(b_1,b_2))  
            # A strand double-crossing a tangle:
            dc_tangle = False
            for tang in tangles:
                if doescross_bool((a_1,inter),tang): # upperhalf
                    if doescross_bool((inter, b_2),tang):
                        dc_tangle = True
                if doescross_bool((b_1,inter),tang): # lower half
                    if doescross_bool((inter,a_2),tang):
                        dc_tangle = True
            # double crossing between strands 
            dc_strand = False
            raw = dic.copy()
            raw.pop(s_pair[0])
            raw.pop(s_pair[1])
            strands = list(raw.values())
            for strand in strands:
                if doescross_bool((a_1,inter),strand): # upperhalf
                    if doescross_bool((inter, b_2),strand):
                        dc_strand = True
                if doescross_bool((b_1,inter),strand): # lower half
                    if doescross_bool((inter,a_2),strand):
                         dc_strand = True
            
            if dc_tangle == False and dc_strand ==False:
                arr.append(s_pair)   
        return arr
            
    def getGenerators(self): #HP
        ''' Returns the list of generators.'''
        return self.p_tangle.getAlgebra(self.is_left, self)
    
    def getIdempotents(self):  # Hsin Pei - come back and edit
        ''' Returns the list of Idempotents, using the corresponding function in `tangle`'''
        return self.p_tangle.getIdempotents(self.is_left)
    
    def _multiplyRaw(self, gen1, gen2):
        ''' If gen1 and gen2 can be multiplied, return the generator that is
        their product. Otherwise, return None.
        gen1 and gen2 are compatible by default.'''
        if self.mod_6(gen1, gen2):
            return None
        else:
            take_right = not self.is_left
            left_strands = gen1.strands.convert_dict(take_right)
            right_strands = gen2.strands.convert_dict(take_right)
            new_pairs = []      
            for k,v in left_strands.items():
                new_pairs.append((k,right_strands[v]))
            new_pairs = tuple(new_pairs)    
            return Simple_Strand(self, self.is_left, new_pairs)
                
    def multiply(self, gen1, gen2):  
        if not isinstance(gen1, Simple_Strand):
            return NotImplemented
        if not isinstance(gen2, Simple_Strand):
            return NotImplemented
        assert gen1.parent == self and gen2.parent == self, \
        "Algebra not compatible."        
        if sorted(gen1.t) != sorted(gen2.s):
            return E0      
        if gen1.sign_seq != gen2.sign_seq:
            return E0
        prod_raw = self._multiplyRaw(gen1, gen2)      
        if prod_raw is None:
            return E0       
        if self.ring is F2:
            return prod_raw.elt()       
        else:
            pass 
    
    def mod_6(self, gen1, gen2): 
        '''returns True if two generators catches mod relations in figure 6 '''
        if self.is_left: # left side of the tangle
            l_strands = gen1.strands.right_converted
            r_strands = gen2.strands.right_converted
        else: # right side of the tangle
            l_strands = gen1.strands.left_converted
            r_strands = gen2.strands.left_converted       
        c_l = gen1.strands.strandCrossing_coord(not self.is_left)
        c_r = gen2.strands.strandCrossing_coord(not self.is_left)   
        l_tangle = [gen1.tangle.orient_right_lhalf, gen1.tangle.orient_left_lhalf]
        r_tangle = [gen2.tangle.orient_right_rhalf, gen2.tangle.orient_left_rhalf]  
#        print("Entering : Mod 6") #cb and remove
#        print("\nl_strands: {0}".format(l_strands))
#        print("\nr_strands: {0}".format(r_strands))
#        print("\nc_l: {0}".format(c_l))
#        print("\nc_r: {0}".format(c_r))
#        print("\nl_tangle {0}".format(l_tangle))
#        print("\nr_tangle {0}".format(r_tangle))        
        return len(self.mult_two_halfs(l_tangle, l_strands, r_tangle, r_strands)) > 0 \
                or len(self.cross_twice(l_strands, r_strands, c_l, c_r)) > 0 

    ## helper methods for mod 6() 
    def mult_two_halfs(self,left_tangle, left_strands, right_tangle, right_strands):
        ''' Multiplies two half tangles T1 and T2, and mod out first two
        relations in PetKova paper Figure 6. 
        left_tangle : [dict1, dict2] right half of tangle T1, array of dictionary objects of coordinates
                      dict1 : orienting left, dict2: orienting right
        right_tangle: [dict1, dict2] left half of tangle T2, a dictionary object of coordinates
        left_strand: right half of tangle T1 strands, a dictionary object of coordinates
        right_strand: left half of tangle T2 strands, a dictionary object of coordinates
        
        pairs is an array object, which elements of form [start tangle, start strand]
        used for algebra multiplication and for other things later( d+, dm)'''     
        pairs = []
        # tangle orienting right      
        t1 = left_tangle[0]
        t2 = right_tangle[0]   
        if(self.is_left): # left side of the tangle
            left_strands = dict_shift(left_strands,True)# shift left_strand to the right
        else:
            right_strands = dict_shift(right_strands,False) # shift left_strand to the right      
    
        for k,v in t1.items():
            for a,b in left_strands.items():
                if doescross((k,v),(a,b)):
                    new_tangle = (v,t2[v]) # on the right
                    new_strand = (b, right_strands[b]) # on the right
                    if doescross(new_tangle, new_strand):
                        pairs.append([k,a])       
        # tangle orienting left
        t1 = left_tangle[1]
        t2 = right_tangle[1]               
        for k,v in t2.items():
            for b,c in right_strands.items():
                if doescross((v,k),(b,c)):
                    new_tangle = (t1[v], v) # on the left        
                    # search for a strand in left_strand, that ends in b
                    for key,value in left_strands.items():
                        if value == b:
                            new_strand = (key,b) # get the left coordinate     
                    if doescross(new_tangle, new_strand):
                        pairs.append([k, c])
        return pairs

    def cross_twice(self,left_strands, right_strands, c_l, c_r): # cb and change it to true or false later if uncessary
        '''Given left_strand, and right_strand, returns double crossings,
        assumes the strands are multiplicable. 
        c_1 is strandCrossing_coord(False) called from the left_strands
        c_2 is the strandCrossing_coord(True) called from right_strands'''
        # check if compatible cb and remove
        if(self.is_left): # left side of  the tangle
            left_strands = dict_shift(left_strands,True)# shift left_strand to the right
            c_l = dict_shift_double(c_l ,True)
        else:
            right_strands = dict_shift(right_strands,False) # shift left_strand to the right 
            c_r = dict_shift_double(c_r ,False)            
        if not (get_range_dict(left_strands) == get_domain_dict(right_strands)):
            raise TypeError ("strand boundaries don't line up")                    
        double_crossings = []
        for k in list(c_l.items()):
            for v in list(c_r.items()):
                if get_end_dict(k) == get_start_dict(v):
                    double_crossings.append(k)
        return double_crossings

class StrandAlgebraElement(Element):
    ''' An element of strand algebra '''    
    def isIdempotent(self): 
        '''Tests whether this element is an idempotent.'''
        for sa, coeff in self.items():
            if not sa.isIdempotent(): # F2 so no need to check coefficients
                return False
        return True
    
    def invertible(self): 
        '''Tests whether this element is an idempotent.'''
        for sd, coeff in self.items():
            if not sd.isIdempotent():
                return False
        return True
    
    def inverse(self):
        ''' Returns the inverse of this element, if it is invertible.'''
        return self
    
############################### TEST CODE ###################################

dic = {(4,4):(3.5,4), (3.5,4):(3,4), (3,2):(3.5,2), (3.5,2):(4,3), (4,2):(3.5,3), (3.5,3):(3,3), (3,1):(3.5,1), (3.5,1):(4,1)}
tang = TANGLE(dic)
strands = (((3,4),(1,2),(2,1)),((3,2),(0,1)))
CT = TangleModule(F2, tang) #CT(Ti)
sd = StrandDiagram(CT,strands)#Strand Diagram 
left_alg = StrandAlgebra(F2, tang, True)
right_alg = StrandAlgebra(F2, tang, False)
a2 = Simple_Strand(right_alg, False, ((1,2),(3,2),(0,4)))
        
#dic = {(4,4):(3.5,4),(3.5,4):(3,4),(4,2):(3.5,3), (3.5,3):(3,3), (3,2):(3,1)}
#tang = TANGLE(dic)
#print(tang.caps)
#print(tang.orient_left_lhalf)
#print(tang.orient_left_rhalf)
#print(tang.orient_right_rhalf)
#print(tang.orient_right_lhalf)
