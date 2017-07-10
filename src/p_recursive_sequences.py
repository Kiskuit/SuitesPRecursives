#!/usr/bin/python2
# -*- coding: utf8 -*-


# ---------------------------------------------------------
# Author : Mathis Caristan (mathis.caristan@gmail.com)
# Date : 19/06/2017
# Class : PRecSequence -- Class to manipulate and compute
#   with p-recursive sequences in Sage/Python
# Copyright 2017 Mathis Caristan
# ---------------------------------------------------------

# General imports
from __future__ import print_function

# from sage.all import *
# from ore_algebra import *
from sage.structure.element import RingElement
from sage.structure.sequence import Sequence
from sage.structure.sage_object import op_EQ,op_NE

# TODO      * un constructeur (fonction séparée) qui fabrique une suite à partir d'une
#           expression sage du genre factorial(n)*2^n + n,
#           * un moyen de calculer des suites du style u(3*n+2) à partir de u(n)...


class PRecursiveSequence(RingElement):
    r"""
    Class to represent P recursive sequences.
    The representation is defined by two parameters, the initial conditions, and recurrence operator (or annihilator).

    EXAMPLES::
        sage: from parent_p_recursive_sequences import *
        sage: P = ParentPRecursiveSequences(ZZ['n'])
        sage: n = P.base_ring().gen()
        sage: Sn = P.generator()
        sage: fibo = PRecSequence (P, [0,1], Sn**2-Sn-1)
        sage: fibo[8] == sloane.A000045[8], fibo[8]
        (True, 21)
        sage: fibo[5:9]
        [5, 8, 13, 21]

    """

    _binary_splitting_threshold = 100

    def __init__ (self, parent, condInit, annihilator):
        r"""
        Initializes a sequence. The initial conditions can be a list or a dictionary.
        If it is a list, it is assumed that the indices start at 0, and are consecutive.
        It requires that the order first conditions are consecutive.

        EXAMPLES::
            sage: from parent_p_recursive_sequences import *
            sage: P = ParentPRecursiveSequences(ZZ['n'])
            sage: n = P.base_ring().gen()
            sage: Sn = P.generator()
            sage: fibo1 = PRecSequence(P, [0,1], Sn**2 - Sn - 1)
            sage: fibo2 = PRecSequence(P, {0:0,1:1}, Sn**2 - Sn -1)
            sage: fibo1 == fibo2
            True
        """

        if isinstance(condInit, dict):
            self.cond = condInit.copy()
        elif isinstance(condInit, list):
            self.cond = {i:val for i,val in enumerate(condInit)}
        else :
            raise ValueError ("condInit must be a list or a dict.")
###        self._i = sorted(self.cond.keys())[0]-1
        
        self.annihilator = parent.ore_algebra().coerce(annihilator)
        if self.annihilator.parent() is not parent.ore_algebra():
            raise ValueError("`annihilator` must be in {}.".format(parent.ore_algebra()))
        # Must have enough initial conditions
        if len (self.cond) < self.annihilator.order() : 
            raise ValueError ("Not enough initial conditions.")
        # order first initial conditions must be consecutive
        elif self.annihilator.order() != 0:
            keys = sorted(self.cond)
            if keys[self.annihilator.order()-1] != keys[0] + self.annihilator.order()-1:
                raise ValueError("{} first conditions must be consecutive.".format(self.annihilator.order()))

        RingElement.__init__(self, parent)

    ###############################################################

    def order (self):
        r"""
        Returns the order of a sequence. The order of a sequence, is the order of its annihilator

        EXAMPLES::
            sage: from parent_p_recursive_sequences import *
            sage: Seqs = ParentPRecursiveSequences(ZZ['n'])
            sage: n = Seqs.generator()
            sage: Sn = P.generator()
            sage: fibo = P([0,1], Sn^2-Sn-1)
            sage: fibo.order()
            2
            sage: u = P([0,2,5,-3], (3*n^2-4)*Sn^4 - 3*Sn^3 + (n+1)*Sn^2 - 5
            sage: u.order()
            4
        """
        return self.annihilator.order()

    ###############################################################

    def __iter__(self):
        r"""
        Iterates infinitely over the sequences values (or until a degenerate values)

        INPUT
        """
        def iterfct():
            i = min(self.cond)
            while True:
                ret = self[i]
                i += 1
                yield ret
        return iterfct()

    ###############################################################

    def list(self, *args):
        r"""
        Returns a list of the values taken by the sequence, in the range chosen by the user.
        
        INPUT::
        - ``args`` (len(args)==1) -- args[0] is the stoping index of the list (the list will start at its first element).
                   (len(args)==2) -- args[0] is the starting index of the list, and args[1] is the stoping index.
                   In both cases, it behaves like the 'range()' function, and last index will not be included.
        
        EXAMPLES::
            sage: import parentSeq
            sage: Seqs = ParentPRecursiveSequences(QQ['n'])
            sage: Sn = P.generator()
            sage: fibo = P([0,1], Sn^2-Sn-1)
            sage: fibo.list(7)
            [0, 1, ..., 8]
            sage: fibo.list(5,9)
            [5, 8, 13, 21]

        TODO::
        - Add opt arguments to chose algo and more
        """
        # TODO provide support for step?
        lowest = min(self.cond.keys())
        if len(args)==0:
            raise TypeError('list() takes at least 2 arguments (1 given)')
        if len(args)==1:
            start = lowest
            stop = args[0]
        elif len(args)==2:
            start = args[0]
            stop = args[1]
        else:
            print ("list() :call to list, too many args :", args)
        # start/stop cannot be lower than lowest index
        if stop < lowest or start < lowest :
            err_str = "Index out of bond, indices cannot be lower than "
            err_str += str(lowest) + "."
            raise IndexError(err_str)
        return self[start:stop]

    ###############################################################

    def __getitem__(self,slice_):
        r"""
        Returns an (list of) element(s) of the sequence.

        INPUT::
        -``slice_`` if this parameter is an integer, it returns the element at index slice_ in the sequence.
                    if it is a Python slice, it returns a list of elements in the slice span

        EXAMPLES::
            TODO

        TODO::
        - Add support for step
        - if cond_init return straight away
        - change or as it does not work well with integer (0 or 1 --> 1)
        - return list or element
        - check if start/stop < lowest index

        """
        # Get start, stop and step params
        try :
            start = slice_.start or min(self.cond.keys())
            step = slice_.step or 1
            stop = slice_.stop
        except AttributeError:
            start = slice_ 
            stop = slice_+1
            step = 1

        if stop < start :
            raise IndexError("Upper index must not be smaller than the lower index")

        if start in self.cond.keys() and start == stop-1:
            return self.cond[start]

        ret = self._condTakenIntoAccount(start, stop, step)
        if len(ret) == 1:
            return ret[0]
        else:
            return ret

    ###############################################################

    def _condTakenIntoAccount (self, start, stop, step, algo='auto'):
        r"""
        Computes and returns the elements of the sequenece.

        INPUT::
        -``start`` starting index of the computation
        -``stop`` stoping index of the computation
        -``step`` only the first out of step elements is returned
        -``algo`` the algorithm used to compute the values, if it is 'auto' (default),
        the function automatically determines which of the two following algo it uses.
        If it is 'to_list', it will reccursively coputes all elements.
        If it is 'bsplit', it will use a binary splitting method
        This option is mostly used for debugging purpose, as it is usually better to let
        the function decide which algo should be used.
        """
        # Setup
        ini = self.cond.copy()
        ord_ = self.order() or 1 # To avoid problems with 0-order sequences
        if start not in ini:
            ini[start] = None
        ini[stop] = None
        cond_keys = sorted(ini)
        init_keys, extra_keys = cond_keys[:ord_], cond_keys[ord_:]
        cond = {k:ini[k] for k in init_keys}
        if stop in init_keys:
            # Enforce type of values via sequences
            try:
                return Sequence([ini[k] for k in range(start,stop)], universe=self.parent().values_ring(), use_sage_types=True)
            except TypeError:
                raise TypeError("Values are not in {}".format(str(self.parent().values_ring())))

        # Computation
        prev = min(start, max(init_keys))
        ret = Sequence([], universe=self.parent().values_ring(), use_sage_types=True)
        print ('condTakenIntoacc() : extra :', extra_keys)
        while extra_keys != [] and extra_keys[0] <= stop:
            key = extra_keys.pop(0)
            if key > start:
                # TODO handle case where a key is in the order terms after start (using min(start,order) proly)
                sta = prev
                algo = 'to_list'
            else:
                sta = key-ord_ + 1
            print('condTaken() : ', cond,key, key+1)
            cond_vals = self._computeElements (cond, sta, key+1, algorithm=algo)
            # If it exists, replace values by user's initial conditions
            if ini[key] is not None:
                cond_vals[-1] = ini[key]
            cond = {i : cond_vals[i-sta] for i in range(sta,key+1)}
            prev = key
            if key > start:
                for elem in cond_vals[:-1]:
                    # Enforce type of values via sequences
                    try :
                        ret.append(elem)
                    except TypeError:
                        raise TypeError("Values must be in {}".format(str(self.parent().values_ring())))
        return ret 

    ###############################################################

    def _computeElements (self, cond, start, stop, algorithm='auto'):
        r"""
        TODO
        """
        min_ = min(cond.keys())
        if algorithm=='auto':
            if start < min_ + self._binary_splitting_threshold:
                algorithm = 'to_list'
            else:
                algorithm = 'bsplit'

        if algorithm == 'to_list':
            cond_vals = [cond[k] for k in sorted(cond.keys())]
        else:
            P,Q = self.annihilator.forward_matrix_bsplit(start-min_,start=min_)
            if Q==0:
                raise Exception ("degenerate values in the sequence.")
            # Why did we use that before, and not anymore..?
            cond_vals = (P*vector([cond[k] for k in sorted(cond.keys())])/Q)
            min_ = start

        ret = self.annihilator.to_list (cond_vals, stop-min_, start=min_)[start-min_:stop-min_]
        return ret

    ###############################################################


    def _add_ (self, other):
        # TODO add cond init in case of discrepancy
        _class = self.__class__
        # Compute new annihilator
        sum_annihilator = self.annihilator.lclm(other.annihilator)
        ord_ = sum_annihilator.order()
        # Set of keys that need to be computed
        key_set = set()
        start = max(self.cond.keys()[0], other.cond.keys()[0])
        for e in range(start, start + ord_): # base init cond
            key_set.add(e)
        for e in self.cond: # additional init cond
            if e > start:
                key_set.add(e)
        leadPol = sum_annihilator[ord_]
        roots = leadPol.roots(multiplicities=False)
        for r in roots: # degenerate vals
            r += ord_
            if r in ZZ and r > start:
                key_set.add(r)
        sum_cond = {}
        for e in key_set:
            try:
                left = self[e]
                right = other[e]
                sum_cond[e] = left+right
            except: # TODO handle exception (which type)
                raise NotImplementedError

        return _class(self.parent(), sum_cond, sum_annihilator)

    ###############################################################

    def _sub_ (self, other):
        # TODO add cond init in case of discrepancy

        _class = self.__class__
        # Compute new annihilator
        sub_annihilator = self.annihilator.lclm(other.annihilator)
        # Compute initial conditions
        min_ = max(self.cond.keys()[0], other.cond.keys()[0])
        sub_cond = {}
        for i in range(min_, min_ + sub_annihilator.order()):
            sub_cond[i] = self[i] - other[i]
        # Compute additional terms to avoid degenerate values
        leadPol = sub_annihilator[sub_annihilator.order()]
        roots = leadPol.roots()
        for r,_ in roots :
            r += sub_annihilator.order()
            if (r >= min_):
                try:
                    left = self[r]
                    right = other[r]
                    sub_cond[r] = left-right
                except:  # TODO handle exception (and which one precisely) better
                    break

        return _class(self.parent(), sub_cond, sub_annihilator)
            
    ###############################################################

    def _mul_ (self, other):
        #TODO add cond in case of discrepancy

        _class = self.__class__
        # Compute new annihilator
        prod_annihilator = self.annihilator.symmetric_product(other.annihilator)
        # Compute new init cond
        prod_cond = {}
        _min = max(self.cond.keys()[0], other.cond.keys()[0])
        for i in range(_min, _min+prod_annihilator.order()):
            prod_cond[i] = self[i] * other[i]
        # Compute additional terms to avoid degenerate values
        leadPol = prod_annihilator[prod_annihilator.order()]
        roots = leadPol.roots()
        for r,_ in roots :
            r += prod_annihilator.order()
            if (r >= min_):
                try:
                    left = self[r]
                    right = other[r]
                    prod_cond[r] = left*right
                except:  # TODO handle exception (and which one precisely) better
                    break

        return _class (self.parent(), prod_cond, prod_annihilator)

    ###############################################################
    
    def _richcmp_ (self, other, op):
        if op == op_EQ or op == op_NE :
            try : 
                sub = self - other
            except TypeError :
                return NotImplementedError
            if sub.is_const() and sub[sub.cond.keys()[0]] == 0:
                return op == op_EQ
            else :
                return op == op_NE
        else : 
            raise NotImplementedError

    ###############################################################

    def is_const(self):
        r"""
        Returns True if the sequence is constant, False otherwise

        TODO::
        - add possibility to test only from a given value upwards
        """
        if self.order() == 0:
            for _,v in self.cond.iteritems():
                if v != 0:
                    return False
            return True
        min_ = min(self.cond)
        cst = self - min_
        min_ = min(cst.cond)
        vals = cst[min_,min_+cst.order()+1]
        print ('is_const :', vals)
        for v in vals:
            if v != 0:
                return False
        return True

    ###############################################################

    def _repr_(self):
        _str = "annihilator : "+str(self.annihilator)+"\n"
        _str += "values : "+str(self.list(9))+" ...\n"
        return "P-recursive sequence\n"+ _str

    ###############################################################
    ###############################################################
    ###############################################################

if __name__ == "__main__" :
    from sage.all import *
    from parent_p_recursive_sequences import ParentPRecursiveSequences
    Seqs = ParentPRecursiveSequences(ZZ['n'])
    Sn = Seqs.generator()
    fibo = Seqs([0,1], Sn**2 - Sn - 1)
    print(fibo + Seqs())
###   
###       # algebre d'Ore en les variable Sn et n
###       A,n = ZZ["n"].objgen()
###       R,Sn = OreAlgebra(A,"Sn").objgen()
###   
###       # algebre d'Ore en les variable Sx et x
###       A2,x = ZZ["x"].objgen()
###       R2,Sx = OreAlgebra(A2,"Sx").objgen()
###       
###       ##                                 ##
###       # Suite avec des racines Dégénerés  #
###       ##                                 ##
###       cond = {0:0, 1:1, 2:1, 3:2, 4:3, 5:5, 7:13}
###       # u1 = (n-1)*(n-2)*Sn**3 - (n-1)*(n-2)*3*Sn - (n-1)*(n-2)*8
###       u1 = (n-1)*(n-2)*(n-203)*Sn**2 - (n-1)*(n-2)*3*Sn - (n-1)*(n-2)*8
###       s1 = PRecSequence (cond, u1)
###       print("s1:")
###       print (s1[20:23])
###       print (s1[103:110])
###       ## fin test suite avec racine dégénéré ##########""
###   
###   
###       ##  fibonacci  ##
###       cond2 = {0:0, 1:1, 2:1, 3:2, 4:3, 5:5}
###       u2 = Sn**2 - Sn - 1
###       fib = PRecSequence (cond2, u2)
###       print("fib[0:10]:")
###       print(fib[0:10])
###   
###       ##  factoriel  ##
###       cond3 = [1]
###       u3 = Sx - x -1
###       fact = PRecSequence (cond3, u3)
###       print("fact[0:10]")
###       print(fact[0:10])
###   
###       ##  somme de fibonacci et factoriel  ##
###       fibfact = fib + fact
###       print("somme (fib+fact)[0:15]")
###       print(fibfact[0:15])
###   
###       ## utilisation de conversion dans mul
###       print("fibonacci * 0.5")
###       print(fib*0.5)
###       ## utilisation de conversion dans add
###       print("somme de fib + 2")
###       print(fib+2)
###   
###       ##  const suites  ##
###       SConst = constPRecSequence(2)
###       Sconst2 = PRecSequence(const = 3)
###   
###       ## ajout de fib et d'une suite constante
###       print("somme de Sconst2 + fib")
###       print(Sconst2 + fib)
###   
###       ##  Test Guessing 1  ##
###       print("test guessing")
###       print(PRecSequence(cond = [1,3,5,7,9,11], Ore = R))         #ok
###       print(PRecSequence(cond = [0,2,4,8,16,32,64], Ore = R))     #ok
###       print(PRecSequence(cond = [0,1,2,3,4,5], Ore = R))          #ok
###       print(PRecSequence(cond = [0,1], Ore = R))                  #fail "no relation found"
###   
