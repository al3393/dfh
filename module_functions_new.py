from tangle_ import Simple_Strand, StrandAlgebra,TANGLE,Strands, Idempotent, StrandDiagram,TangleModule
from utility import mod_helper, mod_between
import numpy as np
import operator
#from algebra import DGAlgebra, Element, Generator, Tensor, TensorGenerator
import statistics as st
from utility import orientation,orientation_i, complement,generate_subset,\
 doescross, doescross_simple, intersections, simple_intersections, F2
from Algebra import DGAlgebra, Element, Generator, Tensor, TensorGenerator
from Algebra import E0
from utility import get_domain, get_range, generate_bijections, combinations, \
                get_domain_dict, get_range_dict, doescross_simple_rc,in_between_list, \
                reorganize_sign, reorganize_sign_2,dict_shift, get_start_dict, get_end_dict,\
                dict_shift_double, doescross_bool, mod_between, mod_helper,F2,conc_strands

def del_l_raw(gen, algebra):
    ''' del_L map that sends CT-(T_i) to  A-(dRT) * I(-dRT) * CT-(T_i)
    gen is an element of CT-(T_i).
    Algebra denotes A-(dRT) here'''
    if not isinstance(gen, StrandDiagram):
        return NotImplemented
    if not isinstance(algebra,StrandAlgebra):
        return NotImplemented
    left_idem = gen.getLeftIdem().toAlgElt(algebra)
    print("+++ IDEMPOTENT!!!! ++++ ")
    print(left_idem)
    print("+++ IDEMPOTENT!!!! ++++ ")
    return d_m_raw_A(left_idem, gen)
    # cb and change
#    if prod_raw == []:
#        return E0
#    if self.ring is F2: # cb and change
#        return prod_raw.elt()   
#    else:
#        pass

def m2(gen, alg_gen):
    ''' m2 map that sends CT-(T) * I(-dRT) * A-(dRT) to CT-(T).
    * gen is a generator in CT-(Ti)
    * alg_gen is a generator of  A-(dRTi) 
    gen1 here refers to the strand diagram generator of right most elementary tangle T_i '''  
    if not isinstance(gen, StrandDiagram):
        return NotImplemented
    if not isinstance(alg_gen, Simple_Strand):
        return NotImplemented
    assert gen.tangle == alg_gen.p_tangle, " Not compatible. Different Tangle."
    assert alg_gen.is_left == False, " Not compatible, Algebra takes -dLT."    
    if not gen.strands.rightCompatible(gen, alg_gen):
        return E0
    prod_raw = m2_raw(gen, alg_gen)
    if prod_raw is None:
        return E0
    if self.ring is F2: # cb and change
        return prod_raw.elt()
    else:
        pass

def m2_raw(gen, alg_gen): # checked
    ''' If gen and alg_gen can be multiplied (i.e.,the strands pairs matching up), return the generator in CT-(T) that is
    their product. Otherwise return None. gen1 and alg_gen 2 are compatible by default.''' 
    if mod_6_m2(gen, alg_gen):
        return None
    else: #concantenate
        l_raw = gen.strands.data[1]
        r_raw = alg_gen.strands.data[0]
        new = conc_strands(l_raw, r_raw)
        return gen.replace_2([[],list(l_raw)],[[],new])
    
def mod_6_m2(gen, alg_gen):
    ''' Mod 6 relations for m2'''    
    l_tangle = [gen.tangle.orient_right_rhalf, gen.tangle.orient_left_rhalf]
    r_tangle = [alg_gen.tangle.orient_right_lhalf, alg_gen.tangle.orient_left_lhalf]  
    l_strands = gen.strands.right_converted
    r_strands = alg_gen.strands.left_converted
    return len(mult_two_halfs(l_tangle, l_strands, r_tangle, r_strands)) > 0
    
def d_plus_raw(gen): #checked
    ''' d_+ map that smooths a black-black crossing of a generator `gen`
    of CT-(T_i)) contained in right halfs ( i- 1/2, i)
    Returns a list of elements of the form ((s1, s2), d_plus_raw_term), where
    s1 < s2 are pairs of strands that d_plus is applied, and
    diff_term is a generator in d_plus() obtained by uncrossing these two
    strands. Together they specify all terms in gen.d_plus(). '''    
    # applying it to line Bj
    lst = []
    l_strands = list(gen.strands.left_converted.items())
    r_strands = list(gen.strands.right_converted.items())
    types = mod_helper(l_strands, r_strands)
    t_r = list(gen.tangle.r_pairs_wc.items())
    r_crossings = types[4] # type 4 : right crossings
    for s1, s2 in r_crossings: 
        r_s_pairs = r_strands.copy()
        r_s_pairs.remove(s1)
        r_s_pairs.remove(s1)
        is_crossed = False
        for strands in r_s_pairs:# check double black - black crossing
            if mod_between(strands, (s1[0],s2[0]),(s1[1],s2[1]), False) == 0:
                is_crossed = True
                break
        for tangle in t_r:# tangle double cross a strand
            if mod_between(tangle, (s1[0],s2[0]),(s1[1],s2[1]), False) == 0:
                is_crossed = True
                break
        if not is_crossed:
                lst.append(mod_helper_2(gen, s1, s2, is_crossed))
    return lst

def d_minus_raw(gen): # checked
    '''d_- map that introduces a black-black crossing to a generator 'gen`
    of CT-(T_i) contained in left halfs, (i, i+1/2).
    Returns a list of elements of the form ((s1, s2), d_minus_raw_term), where
    s1 < s2 are pairs of strands that d_minus is applied, and
        diff_term is a generator in d_minus() obtained by crossing these two
        strands. Together they specify all terms in gen.d_minus(). ''' 
    lst = [] # applying it to line Bj
    l_strands = list(gen.strands.left_converted.items())
    r_strands = list(gen.strands.right_converted.items())
    types = mod_helper(l_strands, r_strands)
    t_l = list(gen.tangle.l_pairs_wc.items())
    # type 1 : left horizontals
    horizontals = types[1]
    for s1, s2 in horizontals:
        l_s_pairs = l_strands.copy() 
        is_crossed = False
        l_s_pairs.remove(s1)
        l_s_pairs.remove(s2)
        for strands in l_s_pairs: # check double strand-strand crossing
            if mod_between(strands, (s1[0],s2[0]),(s1[1],s2[1]),True)==0:
                is_crossed = True  
                break
        for tangle in t_l:# tangle double cross a strand
            if mod_between(tangle, (s1[0],s2[0]),(s1[1],s2[1]),True) == 0:
                is_crossed = True
                break
        if not is_crossed:
                lst.append(mod_helper_2(gen, s1, s2))
    return lst

def mod_helper_3(left_idem, gen, s1, s2): #CB
    ''' Used for `d_m_raw_A`.  Returns the new algebra element (class Simple_Strand)
    and generator(class: Strand Diagram) with index s1 and s2 swapped.
    returns [new_algebra_element, generator]'''

    if isinstance(s1[0][0], int): # i.e a generator element
        sd_1 = gen.strands.get_strand_index(s1, True)
        sd1_left = False
    else: # i.e an algebra element
        sd_1 = left_idem.strands.get_strand_index(s1, False)
        sd1_left = True      
    if isinstance(s2[0][0], int): # i.e a generator element
        sd_2 =  gen.strands.get_strand_index(s2, True)
        sd2_left = False
    else: # i.e an algebra element
        sd_2 = left_idem.strands.get_strand_index(s2, False)
        sd2_left = True
    
    if sd1_left and sd2_left: # two algebra elements
        new_strand = [(sd_1[0],sd_2[1]), (sd_2[0], sd_1[1])]
        return ((sd_1[1],sd_2[1]), left_idem.replace([sd_1, sd_2],new_strand), gen)   
    elif not sd1_left and not sd2_left: # two strand diagram elements
        new_strand = [(sd_1[0],sd_2[1]), (sd_2[0], sd_1[1])]
        return ((sd_1[0],sd_2[0]), left_idem, gen.replace([sd_1,sd_2], new_strand, True))  
    elif sd1_left and not sd2_left: # left algebra element, right strand diagram element
        new_strand_alg = [(sd_1[0], sd_2[0])]
        new_strand_gen = [(sd_1[1],sd_2[1])]
        return ((sd_1[1],sd_2[0]),left_idem.replace([sd_1],new_strand_alg), \
                    gen.replace([sd_2], new_strand_gen, True))
    else: # never occurs # CB and remove
        return NotImplemented

def d_m_raw_A(left_idem, gen):
    '''dm map that picks two pair of points along  A_j, between x and algebra element
    E_L_D(x), a left idempotent of generator x.
    Exchanges two ends of the corresponding pair of black strands.
    if `alg_left` == True, A(-dLT_i) * CT(T_ii)
    Returns a list of elements of the form ((s1, s2), new_alg_element, new_gen), where
    s1 < s2 are pairs of strands that d_m is applied, and
        diff_term is a generator in d_plus() obtained by uncrossing or  crossing
        these two strands. Together they specify all terms in d_m()'''
     # So far only implemented for A(-dLT_i) * CT(T_ii) 
    if not isinstance(left_idem, Simple_Strand):
        raise TypeError("The left element, must be an idempotent in the form of Strand Algebra Element")
    lst =  []  # CB
    l_strands = list(left_idem.strands.right_converted.items())
    r_strands = list(gen.strands.left_converted.items())
    t_l = list(left_idem.tangle.r_pairs_wc.items())
    t_r = list(gen.tangle.l_pairs_wc.items())
    types = mod_helper(l_strands, r_strands)
    # left no cross strands 
    l_horizontals = types[1]
    for s1, s2 in l_horizontals:
#        print("\nCurrent left horiontals: {0} and {1}".format(s1,s2))
        is_crossed = False
        for tangle in t_l: # left half tangles
            if abs(mod_between(tangle, (s1[0],s2[0]),(s1[1],s2[1]), True)) == 1:
                is_crossed = True
                print("tangle left Modded")
                break
        for tangle in t_r: # right half tangles
            x = (tuple(np.add(s1[0],(1,0))), tuple(np.add(s2[0],(1,0))))
            if mod_between(tangle, (s1[1], s2[1]),x, False ) != -2:
                is_crossed = True
                print("tangle right Modded")
                break
        if not is_crossed:
            lst.append(mod_helper_3(left_idem, gen, s1, s2))
    # right cross strands
    r_cross = types[4]
    for s1, s2 in r_cross:
#        print("\nCurrent right cross strands: {0} and {1}".format(s1,s2))
        is_crossed = False
        for tangle in t_l: # left half tangles
            x = (tuple(np.subtract(s1[1],(1,0))), tuple(np.subtract(s2[1],(1,0))))
            if mod_between(tangle,x,(s1[0],s2[0]),True) != -2:
                is_crossed = True
                break
        for tangle in t_r: # right half tangles
            if abs(mod_between(tangle, (s1[0],s2[0]),(s1[1],s2[1]),False)) == 1:
                is_crossed = True
                break
        if not is_crossed:
            lst.append(mod_helper_3(left_idem, gen, s1, s2))
    # left below right
    l_below_r = types[5]
    for s1, s2 in l_below_r:
#        print("\n left below right strands: {0} and {1}".format(s1,s2))
        is_crossed = False
        left_pair = (tuple(np.subtract(s2[0],(0.5,0))), s1[0])
        mid_pair = (s2[0],s1[1])
        right_pair = (tuple(np.add(s1[1],(0.5,0))),s2[1])
        for tangle in t_l: #left half tangles
            if mod_between(tangle, left_pair, mid_pair, True) == -1:
                is_crossed = True
                break
        for tangle in t_r: # right half tangles
            if -2 < mod_between(tangle, mid_pair, right_pair, False) < 1:
                is_crossed = True
                break    
        if not is_crossed:
            lst.append(mod_helper_3(left_idem, gen, s1, s2))
    # left above right strands
    l_above_r = types[6]
    for s1, s2 in l_above_r:
#        print("\n left above right strands: {0} and {1}".format(s1,s2))
        is_crossed = False
        left_pair = (tuple(np.subtract(s2[0],(0.5,0))), s1[0])
        mid_pair = (s2[0],s1[1])
        right_pair = (tuple(np.add(s1[1],(0.5,0))),s2[1])
        for tangle in t_l: # left half tangle
            if mod_between(tangle, left_pair, mid_pair, True) == 1:
                is_crossed = True
                break
        for tangle in t_r: # right half tangle
            if mod_between(tangle, mid_pair, right_pair, False) > -1:
                is_crossed = True
                break
        if not is_crossed:
            lst.append(mod_helper_3(left_idem, gen, s1, s2))  
    return lst

def d_m_raw_B(gen): # checked
    ''' dm map that picks two pair of points along Bj.
    Exchanges two ends of the corresponding pair of black strands of generator
    `gen` in CT(Ti)
        Returns a list of elements of the form ((s1, s2), d_m_raw_term), where
    s1 < s2 are pairs of strands that d_m is applied, and
        diff_term is a generator in d_plus() obtained by uncrossing or  crossing
        these two strands. Together they specify all terms in d_m(). '''
    lst = []
    l_strands = list(gen.strands.left_converted.items())
    r_strands = list(gen.strands.right_converted.items())
    types = mod_helper(l_strands, r_strands)
    t_l = list(gen.tangle.l_pairs_wc.items())
    t_r = list(gen.tangle.r_pairs_wc.items())
    # right no cross strands
    r_horizontals = types[3]
    for s1,s2 in r_horizontals:
        is_crossed = False
        for tangle in t_r: # right half tangles
            if abs(mod_between(tangle, (s1[0],s2[0]),(s1[1],s2[1]),False)) == 1:
                is_crossed = True
                break
        for tangle in t_l: #left half tangles
            x = (tuple(np.subtract(s1[1],(1,0))), tuple(np.subtract(s2[1],(1,0))))
            if mod_between(tangle, x,(s1[0],s2[0]), True) != -2:
                is_crossed = True
                break
        if not is_crossed:
            lst.append(mod_helper_2(gen, s1, s2))
    # left cross strands
    l_cross = types[2]
    for s1,s2 in l_cross:
        is_crossed = False
        for tangle in t_l: #left half tangles
            if abs(mod_between(tangle, (s1[0],s2[0]),(s1[1],s2[1]),True)) == 1:
                is_crossed = True
                break
        for tangle in t_r: # right half tangles
            x = (tuple(np.add(s1[0],(1,0))), tuple(np.add(s2[0],(1,0))))
            if mod_between(tangle,(s1[1],s2[1]),x, False) != -2:
                is_crossed = True 
                break
        if not is_crossed:
            lst.append(mod_helper_2(gen, s1, s2))
    # left above right strands
    l_above_r = types[6]
    for s1, s2 in l_above_r:
        is_crossed = False
        left_pair = (tuple(np.subtract(s2[0],(0.5,0))), s1[0])
        mid_pair = (s2[0],s1[1])
        right_pair = (tuple(np.add(s1[1],(0.5,0))),s2[1])
        for tangle in t_l: # left half tangles
            if -2 < mod_between(tangle, left_pair, mid_pair, True) < 1:
                is_crossed = True
                break
        for tangle in t_r: # right_half tangles
            if mod_between(tangle, mid_pair, right_pair, False) == -1:
                is_crossed = True
                break
        if not is_crossed:
            lst.append(mod_helper_2(gen, s1, s2))
    # left below  right strands
    l_below_r = types[5]      
    for s1,s2 in l_below_r: #s1 on the left, s2 on the right
        is_crossed = False
        left_pair = (tuple(np.subtract(s2[0],(0.5,0))), s1[0])
        mid_pair = (s2[0],s1[1])
        right_pair = (tuple(np.add(s1[1],(0.5,0))),s2[1])
        for tangle in t_l: # left half tangles
            if mod_between(tangle, left_pair, mid_pair, True) > -1 :
                is_crossed = True
                break
        for tangle in t_r:
            if mod_between(tangle, mid_pair, right_pair, False) == 1:
                is_crossed = True
                break
        if not is_crossed:
            lst.append(mod_helper_2(gen, s1, s2))
    return lst

def mod_helper_2(gen,s1, s2):
    '''returns the new generator, used for dm(along Bj, d+, d-  '''
    if isinstance(s1[0][0], int):
        sd_1 = gen.strands.get_strand_index(s1, True)
        sd1_left = True
    else:
        sd_1 = gen.strands.get_strand_index(s1, False)
        sd1_left = False  
    if isinstance(s2[0][0], int):
        sd_2 = gen.strands.get_strand_index(s2, True)
        sd2_left = True
    else:
        sd_2 = gen.strands.get_strand_index(s2, False)
        sd2_left = False    
    
    if sd1_left and sd2_left:
        new_strand = ((sd_1[0],sd_2[1]),(sd_2[0],sd_1[1]))
        return (sd_1[1], sd_2[1]), gen.replace((sd_1, sd_2), new_strand, True)
    elif not sd1_left and not sd2_left:
        new_strand = ((sd_1[0],sd_2[1]),(sd_2[0],sd_1[1]))
        return ((sd_1[0], sd_2[0]), gen.replace((sd_1, sd_2), new_strand, False))
    elif sd1_left and not sd2_left:
        new_strand = [[(sd_1[0],sd_2[0])],[(sd_1[1],sd_2[1])]]
        return ((sd_1[1], sd_2[0]), gen.replace_2([[sd_1],[sd_2]], new_strand))
    else: # almost never occurs  # CB and remove
        return NotImplemented
    
def helper(tang_left, tang_right, pair_1, pair_2, pair_3, is_left, option):
    ''' Helper methods for diff methods. Assumes compatibility. Given a generator x, taking the is_left side, 
    exchanges pairs = (a,b), it seems whether the 
    whether it violates `option`. Refer to options.jpg
    for options.  
    * tang_left is dictionary object of pairs 
    * tang_right is dictionary object of pairs 
    * Pair_1 is the coordinate of two pairs of strand on the left
    * Pair_2 is the coordinate of two pairs of strand in the middle
    * Pair_3 is the coordinate of two pairs of strand in the middle.'''
    mod = False
    if is_left: # only check mod relations 1,2,3
        for k,v in tang_left.items():
            if option ==1:
                if mod_helper((k,v), pair_1, pair_2, is_left) == 1: 
                    mod = True
            elif option == 2:
                if mod_helper((k,v), pair_1, pair_2, is_left) == 0:
                    mod = True
            elif option == 3:
                if mod_helper((k,v), pair_1, pair_2, is_left) == -1: 
                    mod = True
            else:
                raise AssertionError("Something is wrong -- options arent either 1, 2 or 3 ")      
    else: # only check mod relations 4,5,6
        for k,v in tang_right.items():
            if option == 4:
                if mod_helper((k,v), pair_2, pair_3, is_left) == 1:
                    mod = True
            elif option ==5:
                if mod_helper((k,v), pair_2, pair_3, is_left) == 0: 
                    mod = True
            elif option == 6:
                if mod_helper((k,v), pair_2, pair_3, is_left) == -1: 
                    mod = True    
            else:
                raise AssertionError("Something is wrong about `is_left`-- options arent either 4, 5, 6" )
    return mod

## helper methods for mod_6_m2() to apply mod relations in figure 6 PetKova Paper
def mult_two_halfs(left_tangle, left_strands, right_tangle, right_strands):
    ''' Multiplies two half tangles T1 and T2, and mod out first two
    relations in PetKova paper Figure 6. 
    left_tangle : [dict1, dict2] right half of tangle T1, array of dictionary objects of coordinates
                  dict1 : orienting left, dict2: orienting right
    right_tangle: [dict1, dict2] left half of tangle T2, a dictionary object of coordinates
    left_strand: right half of tangle T1 strands, a dictionary object of coordinates
    right_strand: left half of tangle T1 strands, a dictionary object of coordinates
    pairs is an array object, which elements of form [start tangle, start strand]'''
    pairs = []
    # tangle orienting right      
    t1 = left_tangle[0]
    t2 = right_tangle[0]   
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
    
#### TEST CODE ####
#
dic = {(4,4):(3.5,4), (3.5,4):(3,4), (3,2):(3.5,2), (3.5,2):(4,3), (4,2):(3.5,3), (3.5,3):(3,3), (3,1):(3.5,1), (3.5,1):(4,1)}
tang = TANGLE(dic)
strands = (((3,4),(1,2),(2,1)),((3,2),(0,1)))
#strands = (((4,4),(2,2),(3,3)),((1,1),(0,0)))

CT = TangleModule(F2, tang) #CT(Ti)
sd = StrandDiagram(CT,strands) # strand diagram
left_idem = sd.left_idem
print("LEFT IDEM: {0}".format(left_idem))
right_idem = sd.right_idem
print("RIGHT IDEM: {0}".format(right_idem))

left_alg = StrandAlgebra(F2, tang, True)
right_alg = StrandAlgebra(F2, tang, False)
print("D_plus raw: ----------------------------")
d_plus = d_plus_raw(sd)
print(d_minus_raw(sd))
#print("a1 : ----------------------------")
#a1 = Simple_Strand(right_alg, False, ((2,2),(1,4))) # for diagram, look at notebook
#a2 = Simple_Strand(right_alg, False, ((2,3),(1,2)))
#a3 = Simple_Strand(right_alg, False, ((1,0),(2,2)))
#mult_result = m2_raw(sd,a1)
#print(type(mult_result)) # works well!
#print("++++++++++++++++")
#d_m = d_m_raw_B(sd) # works well!
#print("+++++++")
#print(len(d_m[0]))
#print(d_m[0][0])
#print(type(d_m[0][1]))
#print("+++++++")

#print(len(d_m[0][0]))
#print(len(d_m[0][1]))
#print(type(d_m[0][1]))
#print(len(d_m[0])) # works well!
print("++++++++++++++++")
#print("+++++++++TESTING Del_L++++++++")
del_l = del_l_raw(sd, right_alg)

print(len(del_l))
print(del_l[0])


#print(len(del_l))
#print(type(del_l[0][1]))
#for new_gen in del_l:
#    print("size: {0}".format(len(new_gen)))
#    print("pairs:{0}".format(new_gen[0]))
#    print("Type:{0}".format(type(new_gen[1])))
#    print("algebra:{0}".format(new_gen[1]))
#    print("generator:{0}".format(new_gen[2]))
#new_alg_elt = del_l[0][1]
#print(new_alg_elt)
#print(new_alg_elt.tangle)
#print(new_alg_elt.pairs)
#print(new_alg_elt)
#for new_gen in d_m:
#    print("size: {0}".format(len(new_gen)))
#    print("pairs:{0}".format(new_gen[0]))
#    print("Type:{0}".format(type(new_gen[1])))
#    print("strand:{0}".format(new_gen[1]))
#print("1:{0}".format(d_m[0]))
#print("2:{0}".format(d_m[1]))
#print("3:{0}".format(d_m[2]))

# figure 9 PetKova Paper


#p1 = ((1,0.5),(1,4.5))
#p2 = ((1.5,1.5), (1.5, 3.5))
#p3 = ((2,1.5),(2,4.5))
#t_left = {(1,7):(1.5,3), (1,2):(1.5,3)}
#t_right = {(1.5,3):(2,5)}
#
#print(helper(t_left,t_right, p1, p2, p3, False, 6 ))
#         