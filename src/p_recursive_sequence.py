#!/usr/bin/python2
# -*- coding: utf8 -*-


# ---------------------------------------------------------
# Author : Mathis Caristan (mathis.caristan@gmail.com)
# Date : 19/06/2017
# Class : PRecursiveSequence -- Class to manipulate and compute
#   with p-recursive sequences in Sage/Python
# Copyright 2017 Mathis Caristan
# ---------------------------------------------------------

# General imports
from __future__ import print_function

# Sage imports
from sage.structure.element import RingElement
from sage.structure.sequence import Sequence
from sage.structure.sage_object import op_EQ,op_NE
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.rings.semirings.non_negative_integer_semiring import NN
from sage.rings.infinity import Infinity
from sage.modules.free_module_element import vector

# TODO      * un constructeur qui fabrique une suite à partir d'une
#           expression sage du genre factorial(n)*2^n + n,
#           * un moyen de calculer des suites du style u(3*n+2) à partir de u(n)...


class PRecursiveSequence(RingElement):
    r"""
    Class to represent P recursive sequences.
    The representation is defined by two parameters, the initial conditions, and recurrence operator (or annihilator).

    EXAMPLES::

        sage: from parent_p_recursive_sequences import *
        sage: P = PRecursiveSequences(ZZ['n'])
        sage: n = P.base_ring().gen()
        sage: Sn = P.shift_operator()
        sage: fibo = PRecursiveSequence (P, [0,1], Sn**2-Sn-1)
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
            sage: P = PRecursiveSequences(ZZ['n'])
            sage: n = P.base_ring().gen()
            sage: Sn = P.shift_operator()
            sage: fibo1 = PRecursiveSequence(P, [0,1], Sn**2 - Sn - 1)
            sage: fibo2 = PRecursiveSequence(P, {0:0,1:1}, Sn**2 - Sn -1)
            sage: fibo1 == fibo2
            True
        """

        # Setup sage
        RingElement.__init__(self, parent)
        # Parameters type checking
        if isinstance(condInit, dict):
            self.cond = condInit.copy()
        elif isinstance(condInit, list):
            self.cond = {i:val for i,val in enumerate(condInit)}
        else :
            raise ValueError ("condInit must be a list or a dict.")
        if len(self.cond) == 0:
            self.cond = {0:0}
        # Annihilator setup
        self._annihilator = parent.ore_algebra().coerce(annihilator)
        ord_ = self.order()
        # Initial conditions setup
        if len (self.cond) < ord_ : 
            raise ValueError ("Not enough initial conditions.")
        elif ord_ > 1:
            keys = sorted(self.cond)
            if keys[ord_-1] != keys[0] + ord_-1:
                raise ValueError("{} first conditions must be consecutive.".format(ord_))
        singular = self.singular_indices()
        for e in sorted(self.cond)[(ord_ or 1):]:
            if e not in singular:
                real_val = self[e]
                if real_val != self.cond[e]:
                    raise ValueError("You provided a wrong value for a non singular term.")
                else:
                    del self.cond[e]


    ###############################################################

    def annihilator(self):
        r"""
        Returns the annihilator of the sequence.
        """
        return self._annihilator
    ###############################################################

    def order (self):
        r"""
        Returns the order of a sequence. The order of a sequence, is the order of its annihilator.

        EXAMPLES::

            sage: from parent_p_recursive_sequences import *
            sage: Seqs = PRecursiveSequences(ZZ['n'])
            sage: n = Seqs.base_ring().gen()
            sage: Sn = Seqs.shift_operator()
            sage: fibo = Seqs([0,1], Sn^2-Sn-1)
            sage: fibo.order()
            2
            sage: u = Seqs([0,2,5,-3], (3*n^2-4)*Sn^4 - 3*Sn^3 + (n+1)*Sn^2 - 5)
            sage: u.order()
            4
        """
        return self.annihilator().order()

    ###############################################################

    def __iter__(self):
        r"""
        Iterates indefinitely over the sequence's values (or until a degenerate value appears)

        INPUT
        """
        # TODO start at a given point
        def iterfct():
            ord_ = self.order()
            cond = [(k,self.cond[k]) for k in sorted(self.cond)][:ord_]
            for k,v in cond:
                i = k+1
                yield v
            while True:
                if i in self.cond:
                    cond.append((i,self.cond[i]))
                    cond = cond[1:]
                    i += 1
                    yield cond[-1][1]
                else:
                    ret = self._computeElements (dict(cond), i, i+1,algorithm="to_list")[0]
                    if ret is None:
                        raise Exception ("Degenerate value")
                    cond.append((i,ret))
                    cond = cond[1:]
                    i += 1
                    yield ret
        return iterfct()

    ###############################################################

    def list(self, *args):
        r"""
        Returns a list of the values taken by the sequence, in the range chosen by the user.
        
        INPUT:

        - ``args`` (len(args)==1) -- args[0] is the stoping index of the list (the list will start at its first element).
                   (len(args)==2) -- args[0] is the starting index of the list, and args[1] is the stoping index.
                   In both cases, it behaves like the 'range()' function, and last index will not be included.
        
        EXAMPLES::

            sage: from parent_p_recursive_sequences import *
            sage: Seqs = PRecursiveSequences(QQ['n'])
            sage: Sn = Seqs.shift_operator()
            sage: fibo = Seqs([0,1], Sn^2-Sn-1)
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
            start = min(self.cond)
            stop = start + 10
        elif len(args)==1:
            start = lowest
            stop = args[0]
        elif len(args)==2:
            start = args[0]
            stop = args[1]
        else:
            print ("list() :call to list, too many args :", args)
            raise TypeError("Too many arguments for list")
        # start/stop cannot be lower than lowest index
        if stop < lowest or start < lowest :
            err_str = "Index out of bounds, indices cannot be lower than "
            err_str += str(lowest) + "."
            raise IndexError(err_str)
        ret = self[start:stop]
        try:
            return list(ret)
        except TypeError:
            return [ret]

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
        - return list or element
        - check if start/stop < lowest index

        """
        # Get start, stop and step params
        try :
            if slice_.start is None:
                start = min(self.cond)
            else:
                start = slice_.start
            if slice_.step is None:
                step = 1
            else:
                step = slice_.step
            stop = slice_.stop
        except AttributeError:
            start = slice_ 
            stop = slice_+1
            step = 1

        if stop < start :
            raise IndexError("Upper index must not be smaller than the lower index")

        if start in self.singular_indices()  and start == stop-1:
            try:
                return self.cond[start]
            except KeyError:
                raise ValueError("You must provide a value for this degenerate index")

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
        -``stop`` stopping index of the computation
        -``step`` only the first out of step elements is returned
        -``algo`` the algorithm used to compute the values, if it is 'auto' (default),
        the function automatically determines which of the two following algo it uses.
        If it is 'to_list', it will reccursively coputes all elements.
        If it is 'bsplit', it will use a binary splitting method
        This option is mostly used for debugging purposes, as it is usually better to let
        the function decide which algo should be used.
        """
        # Setup
        ini = self.cond.copy()
        singular = self.singular_indices()
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
                return Sequence([ini[k] for k in range(start,stop)],
                             universe=self.parent().values_ring(), use_sage_types=True)
            except TypeError:
                raise TypeError("Values are not in {}".format(str(self.parent().values_ring())))
        # Computation
        prev = min(start, max(init_keys))
        ret = Sequence([], universe=self.parent().values_ring(), use_sage_types=True)
        while extra_keys != [] and extra_keys[0] <= stop:
            key = extra_keys.pop(0)
            if key > start:
                sta = prev
                algo = 'to_list'
            else:
                sta = key-ord_ + 1
            cond_vals = self._computeElements (cond, sta, key+1, algorithm=algo)
            # If it is a singular value, and its value is set
            if key in singular and ini[key] is not None:
                cond_vals[-1] = ini[key]
            # If it is not, check if the value provided by user is correct
            ###   elif key not in singular and key != stop:
            ###       if ini[key] != cond_vals[-1]:
            ###           raise ValueError("You provided a wrong value for a term that is not singular")
            cond = {i : cond_vals[i-sta] for i in range(sta,key+1)}
            prev = key
            if key > start:
                for elem in cond_vals[:-1]:
                    if elem is None :
                        ret = ret + [None]
                    else:
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
            P,Q = self.annihilator().forward_matrix_bsplit(start-min_,start=min_)
            if Q==0:
                raise Exception ("degenerate values in the sequence.")
            # Why did we use that before, and not anymore..?
            cond_vals = (P*vector([cond[k] for k in sorted(cond.keys())])/Q)
            min_ = start

        ret = self.annihilator().to_list (cond_vals, stop-min_, start=min_)[start-min_:stop-min_]
        return ret

    ###############################################################

    def _operation(self, other, op_terms, op_annihil):
        _class = self.__class__
        #new_annihilator = self.annihilator().op_annihil(other.annihilator())
        new_annihilator = op_annihil(other.annihilator())
        new_annihilator = self.parent().ore_algebra()(new_annihilator)
        ord_ = new_annihilator.order()
        key_set = set()
        start = max(self.cond.keys()[0], other.cond.keys()[0])
        # ord_ first keys
        for e  in range(start, start + ord_):
            key_set.add(e)
        leadPol = new_annihilator[ord_]
        roots = leadPol.roots(multiplicities=False)
        # roots
        for r in roots:
            r+= ord_
            if r in ZZ and r > start:
                key_set.add(r)
        new_cond = {}
        for e in key_set:
            try:
                left = self[e]
                right = other[e]
                new_cond[e] = op_terms(left,right)
            except: # TODO handle exception (which type)
                continue
        return _class(self.parent(), new_cond, new_annihilator)

    ###############################################################

    def _add_ (self, other):
        op = self.annihilator().lclm
        return self._operation(other, lambda x,y : x+y, op)

    def add (self, other, method='default'):
        if method == 'default':
            return self + other
        elif method == 'guess':
            pass
        else: # TODO
            raise ValueError("method must either be `default` or `guess`")

    ###############################################################

    def _sub_ (self, other):
        op = self.annihilator().lclm
        return self._operation(other, lambda x,y : x-y, op)
            
    ###############################################################

    def _mul_ (self, other):
        op = self.annihilator().symmetric_product
        return self._operation(other, lambda x,y : x*y, op)

    ###############################################################

    @cached_method
    def __nonzero__(self):
        min_ = min(self.cond)
        if self.is_const() and self.cond[min_]==0:
            return False
        return True

    ###############################################################

    @cached_method
    def singular_indices (self):
        ord_ = self.order()
        pol = self._annihilator[ord_]
        roots = pol.roots(multiplicities=False)
        res = []
        for r in roots:
            if r in ZZ and r+ord_ > min(self.cond):
                res.append(r+ord_)
        return res

    ###############################################################
    
    def _richcmp_ (self, other, op):
        if op == op_EQ or op == op_NE :
            minSelf = min(self.cond)
            minOther = min(other.cond)
            try : 
                sub = self - other
            except TypeError :
                raise NotImplementedError
            if minSelf == minOther and sub[sub.cond.keys()[0]] == 0 and sub.is_const():
                return op == op_EQ
            else :
                return op == op_NE
        else : 
            raise NotImplementedError

    ###############################################################

    @cached_method
    def is_const(self):
        r"""
        Returns True if the sequence is constant, False otherwise

        EXAMPLE::

            sage: import parent_p_recursive_sequences as par
            sage: Seqs = par.PRecursiveSequences(ZZ['n']);
            sage: u,v = Seqs.an_element(), Seqs.zero()
            sage: u.is_const(), v.is_const()
            (False, True)

        TODO::

        - add possibility to test only from a given value upwards
        """
        # If order is 0, the sequence is 0 everywhere but the initial conditions
        if self.order() == 0:
            for _,v in self.cond.iteritems():
                if v != 0:
                    return False
            return True
        min_ = min(self.cond)
        cst = self - self.cond[min_]
        min_ = min(cst.cond)
        vals = cst.list(min_,min_+cst.order())
        for v in vals:
            if v != 0:
                return False
        return True

    ###############################################################

    def _repr_(self):
        try:
            vals = self.list()
        except ValueError:
            return "This sequence contains a singuar value"
        str_ = "["
        for v in vals:
            str_ += str(v) + ", "
        str_ += "...]"
        min_ = min(self.cond.keys())
        if min_ != 0:
            str_ += " (starting at n={})".format(min_)
        return str_

    ###############################################################
    ###############################################################
    ###############################################################

if __name__ == "__main__" :
    #from sage.all import *
    from parent_p_recursive_sequences import PRecursiveSequences
    Seqs = PRecursiveSequences(ZZ['n'])
    Sn = Seqs.shift_operator()
    u = PRecursiveSequence (Seqs, [], 1)
    u.is_zero()
###       cond3 = [1]
###       u3 = Sx - x -1
###       fact = PRecursiveSequence (cond3, u3)
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
###       SConst = constPRecursiveSequence(2)
###       Sconst2 = PRecursiveSequence(const = 3)
###   
###       ## ajout de fib et d'une suite constante
###       print("somme de Sconst2 + fib")
###       print(Sconst2 + fib)
###   
###       ##  Test Guessing 1  ##
###       print("test guessing")
###       print(PRecursiveSequence(cond = [1,3,5,7,9,11], Ore = R))         #ok
###       print(PRecursiveSequence(cond = [0,2,4,8,16,32,64], Ore = R))     #ok
###       print(PRecursiveSequence(cond = [0,1,2,3,4,5], Ore = R))          #ok
###       print(PRecursiveSequence(cond = [0,1], Ore = R))                  #fail "no relation found"
###   
