#!/usr/bin/python2
# -*- coding: utf8 -*-

# ---------------------------------------------------------
# Author : Mathis Caristan (mathis.caristan@gmail.com)
# Date : 19/06/2017
# Class : ParentSeq
# Copyright 2017 Mathis Caristan
# ---------------------------------------------------------

from __future__ import print_function

from sage.rings.ring import Ring
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.rings import Rings
from sage.structure.coerce import py_scalar_to_element
#from sage.all import *
from ore_algebra import OreAlgebra, guess
from seq import PRecSequence

class ParentPRecursiveSequences (Ring, UniqueRepresentation):
    r"""
    Ensemble of P recursive sequences which share the same ring for the coefficient,
    and the same ring for the values.

    EXAMPLES::
        sage: from ore_algebra import *
        sage: Seqs = ParentPRecursiveSequences(ZZ['n'])
        sage: map(lambda x:Seqs.has_coerce_map_from(x), [ZZ,QQ,RR])
        [True, False, False]
        sage: Seqs = ParentPRecursiveSequences(QQ['n'], values_ring=CC)
        sage: map(lambda x:Seqs.has_coerce_map_from(x), [ZZ,QQ,RR,CC])
        [True,True,True,True]
    """
    Element = PRecSequence
    # TODO generator as optional arg and if none : gen = 'S'+base.gen() ?
    def __init__(self, base_ring, generator=None, values_ring=None, category=None):
        r"""
        Initializes "self".
        """
        # Default values_ring
        if values_ring is None:
            values_ring = base_ring.base_ring()
        elif not base_ring.base_ring().is_subring(values_ring):
            raise ValueError("`base_ring` must be a part of `values_ring`")
        self._values_ring = values_ring
        # Default generator
        if generator is None:
            generator = 'S' + str(base_ring.gen())

        # TODO catch error, to raise mine
        self._ore_algebra, self._generator = OreAlgebra (base_ring, generator).objgen()
        Ring.__init__(self, base_ring, category=category or Rings())

    def values_ring(self):
        r"""
        Returns the ring of the values.
        """
        return self._values_ring

    def ore_algebra(self):
        r"""
        Returns the ore algebra in which the annihilator of the sequences are.
        """
        return self._ore_algebra

    def generator(self):
        r"""
        Returns the generator of the ore algebra of the sequences.
        Since we work with P recursive sequences, the generators are Standard shift operators (Sn, Sx...).
        """
        return self._generator

    def base_ring(self):
        r"""
        Returns the (polynomial) ring  of the sequences.
        """
        return self.ore_algebra().base_ring()

    def _repr_(self):
        r"""
        Returns the representation of this object.
        """
        ret = "P-Recursive sequences over `{}`, with values in `{}`".format(self.base_ring(),
                                                                            self.values_ring())
        return ret

    def _element_constructor_(self, *args, **kwargs):
        r"""
        Constructs and returns a P recusive sequences.
        If a single argument is given, it tries to interpret it as a sequence (see EXAMPLES).
        Otherwise, it works as the constructor of the element, without the need to provide the parent.


        EXAMPLES::
            sage: Seqs = ParentPRecursiveSequences(ZZ['n'])
            sage: u = Seqs(1)
            [1,1,1, ..., 1...]
            sage: fibo = Seqs({0:1,1:1}, Sn^2-Sn-1)
            sage: fibo[0:15]
            [0, 1, 1, 2, ..., 233, 377]

        """
        if len(args) != 1:
            return self.element_class (self, *args, **kwargs)
        x = py_scalar_to_element(args[0])
        try :
            P = x.parent()
        except AttributeError :
            return self.element_class (self, x, **kwargs)

        # x is in base_ring()
        # TODO take into account cases where there are roots to the pol
        # TODO lower boundary guessing --> `8` hard coded but ugly
        if x in self.base_ring() :
            gen = self.base_ring().gen()
            condInit = [x.subs({gen:0})]
            annihil = -guess ([x.subs({gen:i}) for i in range(8)], self._ore_algebra)
            return self.element_class (self, condInit, annihil)

        # x is a constant of values_ring
        #   return with annihilator Sn - 1 (where Sn is replaced by actual generator)
        if x in self._values_ring :
            return self.element_class (self, {0:x}, self._ore_algebra.gen() - 1)

        # Default case
        return self.element_class (self, x, **kwargs)

    def _coerce_map_from_ (self, S):
        r"""
        Defines the ensembles sage can coerce elements from.
        Coercions exist anuthing that can be coerced into the base ring of the underlying ore algebra,
        and from anything that can be coerced into the values ring.
        The coercion creates constant sequences.
        """
        if self.ore_algebra().base_ring().has_coerce_map_from (S):
            return True
        if self._values_ring.has_coerce_map_from (S):
            return True
        
        # Default case 
        return False

    def guess (self, vals) :
        r"""
        Tries to guess the recurrence of the given list, in the context of "self.ore_algebra()".
        /!\ It uses the `guess` function from ore_algebra, which is currently broken in some cases.
        
        TODO:
        Allow guessing with elements that are not in the ore algebra.
        For instance, [1/2, 1/2, 1/2, 1/2, 1/2...] should give an answer, if the values ring includes 1/2.
        """
        annihil = -guess (vals, self._ore_algebra)
        return self.element_class (self, vals, annihil)



###   if __name__ == "__main__":
###       R = ZZ['n']
###       par = ParentPRecursiveSequences (R,'Sn',RR)
###       print ("-"*32)
###       print (par)
###       print ("-"*32)
###       try:
###           par = ParentPRecursiveSequences (R,'Sx',RR)
###           print (par)
###       except TypeError:
###           print ("Sx is not a generator of the OreAlgebra over {}.".format(R))
###       print ("-"*32)
###       try:
###           par = ParentPRecursiveSequences (ZZ, 'Sn', RR)
###           print (par)
###       except TypeError:
###           print ("ZZ isnt a proper ring for OreAlgebra.")
###       print ("-"*32)
