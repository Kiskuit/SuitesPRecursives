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

from sage.all import *
from ore_algebra import *
from sage.structure.sage_object import op_EQ,op_NE

# TODO      * un constructeur (fonction séparée) qui fabrique une suite à partir d'une
#           expression sage du genre factorial(n)*2^n + n,
#           * un moyen de calculer des suites du style u(3*n+2) à partir de u(n)...


class PRecSequence(RingElement):
    """
    TODO doc
    """

    THRESHOLD = 100

    def __init__ (self, parent, condInit, annihilator, method='auto'):
        """
        """

        if isinstance(condInit, dict):
            self.cond = condInit.copy()
        elif isinstance(condInit, list):
            self.cond = {i:val for i,val in enumerate(condInit)}
        else :
            raise ValueError ("condInit must be a list or a dict")
        self._i = sorted(self.cond.keys())[0]-1
        
        if annihilator not in parent.ore_algebra():
            return ValueError("`annihilator must be in {}.".format(parent.ore_algebra()))
        self.annihilator = annihilator
        if len (condInit) < annihilator.order() : 
            raise ValueError ("Not enough initial conditions")
        if method not in ('bsplit', 'to_list', 'auto'):
            raise ValueError('Unknown value for `method` parameter')
        self._method = method

        RingElement.__init__(self, parent)

###############################################################

    def order (self):
        return self.annihilator.order()

###############################################################

    def __iter__(self):
        return self

###############################################################

    def next(self):
        self._i += 1
        # TODO check if [0] or [-1]
        return (self[self._i])

###############################################################

    def to_list(self, stop, start=None):
        """
        TODO doc
        """
        # TODO provide support for step?

        lowest = min(self.cond.keys())
        if not start :
            start = lowest
        # start/stop cannot be lower than lowest index
        if stop < lowest or start < lowest :
            err_str = "Index out of bond, indices cannot be lower than "
            err_str += str(lowest) + "."
            raise IndexError(err_str)
        return self[start:stop]

###############################################################

    def __getitem__(self,slice_):
        # TODO Add support for step
        # TODO if cond_init return straight away
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

# Abandon of this segment as it is inconsistent in case of degenerated values
# TODO : see if there's a situation where it makes sense
###       def _condJustReplace (self, start, stop, step):
###           # TODO sort things out with extra conds and degenerated vals
###           ret = []
###           cond = [self.cond[i] for i in sorted (self.cond.keys())][:self.order()]
###           if start > self.THRESHOLD + max(self.cond.keys()):
###               # use of forward matrix bsplit to compute ret first terms
###               P,Q = self.annihilator.forward_matrix_bsplit (start)
###               if Q==0:
###                   # TODO find a better exception.
###                   raise Exception ("Degenerated values in the sequence.")
###               # TODO try transpose instead of this blob
###               ret = [e[0] for e in (P*Matrix([[f] for f in cond]))/Q]
###                   
###           else :
###               # Recursive method to compute ret first terms
###               ret = self.annihilator.to_list(cond, start+self.order(), start=min(self.cond.keys()))
###           self.annihilator.to_list(ret, stop-start, start=start, append=True)
###   
###           # Replace with conds
###           for key in self.cond.keys() :
###               if start <= key and key <= stop:
###                   ret[key-start] = self.cond[key]
###   
###           ret = ret[:stop-start:step]
###           if None in ret :
###               raise Exception("Degenerated values in the sequence.")
###           if len(ret) == 1:
###               return ret[0]
###           else :
###               return ret

###############################################################

    def _condTakenIntoAccount (self, start, stop, step):
        cond_keys = [k for k in sorted(self.cond.keys())]
        init_keys, extra_keys = cond_keys[:self.order()], cond_keys[self.order():]
        cond = {k:self.cond[k] for k in init_keys}
        key = None
        # Extra conditions before start
        while extra_keys != [] and extra_keys[0] <= start:
            key = extra_keys.pop(0)
            cond_vals = self._computeElements (cond, key-self.order()+1, key+1, method=self._method)
            cond_vals[-1] = self.cond[key]
            cond = {key-self.order()+1+i:cond_vals[i] for i in range(self.order())}
        # Start to min(extra conditions, stop)
        ret = []
        try:
            key = extra_keys.pop(0)
            assert(key < stop)
            cond_vals = self._computeElements (cond, start, key+1, method = self._method)
            cond_vals[-1] = self.cond[key]
            ret, cond_vals = cond_vals, cond_vals[-self.order():]
            cond = {key-self.order()+1+i:cond_vals[i] for i in range(self.order())}
        except IndexError:
            pass
        except AssertionError:
            ret = self._computeElements (cond, start, stop, method = self._method)
        # Extra conditions after start and before stop
        while extra_keys != [] and extra_keys < stop:
            key = extra_keys.pop(0)
            cond_vals = self._computeElements (cond, key-self.order()+1, key+1, method = self._method)
            cond_vals[-1] = self.cond[key]
            cond = {key-self.order()+1+i:cond_vals[i] for i in range(self.order())}
        if key and start <= key:
            ret += self._computeElements(cond, key+1,stop)
        else:
            ret += self._computeElements(cond, start, stop, method = self._method)
        return ret

###############################################################

    def _computeElements (self, cond, start, stop, method='auto'):
        min_ = min(cond.keys())
        if method=='auto':
            if start < min_+100:
                method = 'to_list'
            else:
                method = 'bsplit'

        if method == 'to_list':
            cond_vals = [cond[k] for k in sorted(cond.keys())]
        else:
            P,Q = self.annihilator.forward_matrix_bsplit(start-min_,start=min_)
            if Q==0:
                raise Exception ("Degenerated values in the sequence.")
            # Why did we use that before, and not anymore..?
            #cond_vals = ((P*Matrix([[f] for f in cond]))/Q).transpose()[0]
            cond_vals = (P/Q).transpose()[-1]
            min_ = start

        ret = self.annihilator.to_list (cond_vals, stop-min_, start=min_)[start-min_:stop-min_]
        return ret

###############################################################


    def _add_ (self, other):
        """
        """
        # TODO add cond init in case of discrepancy

        _class = self.__class__
        # Compute new annihilator
        sum_annihilator = self.annihilator.lclm(other.annihilator)
        # Compute initial conditions
        min_ = max(self.cond.keys()[0], other.cond.keys()[0])
        sum_cond = {}
        for i in range(min_, min_ + sum_annihilator.order()):
            sum_cond[i] = self[i] + other[i]
        # Compute additional terms to avoid degenerated values
        leadPol = sum_annihilator[sum_annihilator.order()]
        roots = leadPol.roots()
        for r,_ in roots :
            r += sum_annihilator.order()
            if (r >= min_):
                try:
                    left = self[r]
                    right = other[r]
                    sum_cond[r] = left+right
                except:  # TODO handle exception (and which one precisely) better
                    break
        return _class(self.parent(), sum_cond, sum_annihilator)
            
###############################################################

    def _sub_ (self, other):
        """
        """
        # TODO add cond init in case of discrepancy

        _class = self.__class__
        # Compute new annihilator
        sub_annihilator = self.annihilator.lclm(other.annihilator)
        # Compute initial conditions
        min_ = max(self.cond.keys()[0], other.cond.keys()[0])
        sub_cond = {}
        for i in range(min_, min_ + sub_annihilator.order()):
            sub_cond[i] = self[i] - other[i]
        # Compute additional terms to avoid degenerated values
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
        """
        """
        #TODO add cond in case of discrepancy

        _class = self.__class__
        # Compute new annihilator
        prod_annihilator = self.annihilator.symmetric_product(other.annihilator)
        # Compute new init cond
        prod_cond = {}
        _min = max(self.cond.keys()[0], other.cond.keys()[0])
        for i in range(_min, _min+prod_annihilator.order()):
            prod_cond[i] = self[i] * other[i]
        # Compute additional terms to avoid degenerated values
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

###############################################################

###       def _eq_(self, other):
###           try :
###               sub = self - other
###               print('hey')
###           except :
###               # No comparison is implemented between the types of self & other
###               print('hey2')
###               return NotImplemented
###   
###           print (type(sub),sub.parent())
###           if sub.is_const() and sub[sub.cond.keys()[0]] == 0:
###               return True
###           return False
###   
###       def _ne_ (self, other):
###           ret = self == other
###           if ret == NotImplemented : 
###               return ret
###           return not self == other

###############################################################

    def is_const(self):
        return True
        for i in self.cond.values():
            if(i != self.cond.values()[0]):
                return False
        #compute enough value
        tab = self.to_list(2*self.order,self.order)
        print(tab)
        for elt in tab:
            if(elt != self.cond.values()[0]):
                return False
        #if is const change self with a reduction????

        #------------------

        return True

###############################################################

    def _repr_(self):
        _str = "annihilator : "+str(self.annihilator)+"\n"
        _str += "values : "+str(self.to_list(9))+" ...\n"
        return "P-recursive sequence\n"+ _str

###############################################################

if __name__ == "__main__" :

    # algebre d'Ore en les variable Sn et n
    A,n = ZZ["n"].objgen()
    R,Sn = OreAlgebra(A,"Sn").objgen()

    # algebre d'Ore en les variable Sx et x
    A2,x = ZZ["x"].objgen()
    R2,Sx = OreAlgebra(A2,"Sx").objgen()
    
    ##                                 ##
    # Suite avec des racines Dégénerés  #
    ##                                 ##
    cond = {0:0, 1:1, 2:1, 3:2, 4:3, 5:5, 7:13}
    # u1 = (n-1)*(n-2)*Sn**3 - (n-1)*(n-2)*3*Sn - (n-1)*(n-2)*8
    u1 = (n-1)*(n-2)*(n-203)*Sn**2 - (n-1)*(n-2)*3*Sn - (n-1)*(n-2)*8
    s1 = PRecSequence (cond, u1)
    print("s1:")
    print (s1[20:23])
    print (s1[103:110])
    ## fin test suite avec racine dégénéré ##########""


    ##  fibonacci  ##
    cond2 = {0:0, 1:1, 2:1, 3:2, 4:3, 5:5}
    u2 = Sn**2 - Sn - 1
    fib = PRecSequence (cond2, u2)
    print("fib[0:10]:")
    print(fib[0:10])

    ##  factoriel  ##
    cond3 = [1]
    u3 = Sx - x -1
    fact = PRecSequence (cond3, u3)
    print("fact[0:10]")
    print(fact[0:10])

    ##  somme de fibonacci et factoriel  ##
    fibfact = fib + fact
    print("somme (fib+fact)[0:15]")
    print(fibfact[0:15])

    ## utilisation de conversion dans mul
    print("fibonacci * 0.5")
    print(fib*0.5)
    ## utilisation de conversion dans add
    print("somme de fib + 2")
    print(fib+2)

    ##  const suites  ##
    SConst = constPRecSequence(2)
    Sconst2 = PRecSequence(const = 3)

    ## ajout de fib et d'une suite constante
    print("somme de Sconst2 + fib")
    print(Sconst2 + fib)

    ##  Test Guessing 1  ##
    print("test guessing")
    print(PRecSequence(cond = [1,3,5,7,9,11], Ore = R))         #ok
    print(PRecSequence(cond = [0,2,4,8,16,32,64], Ore = R))     #ok
    print(PRecSequence(cond = [0,1,2,3,4,5], Ore = R))          #ok
    print(PRecSequence(cond = [0,1], Ore = R))                  #fail "no relation found"

