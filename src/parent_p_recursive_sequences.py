#!/usr/bin/python2
# -*- coding: utf8 -*-

# ---------------------------------------------------------
# Author : Mathis Caristan (mathis.caristan@gmail.com)
# Date : 19/06/2017
# Class : ParentSeq
# Copyright 2017 Mathis Caristan
# ---------------------------------------------------------

from __future__ import print_function

# Sage imports
from sage.rings.ring import Ring
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.rings import Rings
from sage.structure.coerce import py_scalar_to_element
from sage.functions.other import binomial

# Personnal imports
from ore_algebra import OreAlgebra, guess
from p_recursive_sequences import PRecursiveSequence

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
        [True, True, True, True]

    TESTS::
        # Parent creation, unit/zero
        sage: Seqs = ParentPRecursiveSequences(QQ['n']); Sn = Seqs.generator(); n = Seqs.base_ring().gen()
        sage: z1 = Seqs(); z1
        [0, 0, ..., 0, ...]
        sage: z2 = Seqs(0); z2
        [0, 0, ..., 0, ...]
        sage: z1 == z2, z1.is_zero(), z2.is_zero(), z1.is_const(), z2.is_const()
        (True, True, True, True, True)
        sage: Seqs(1), Seqs(1) == Seqs.one(), Seqs.one() == 1
        ([1, 1, ..., 1, ...], True, True)
        sage: Seqs(1).is_one(), Seqs(1).is_const()
        (True, True)
        sage: from random import randint
        sage: annihil = Seqs.ore_algebra().random_element(); ord_ = annihil.order()
        sage: uRnd = Seqs([randint(-20,20) for _ in range(ord_)], annihil)
        sage: (uRnd - uRnd).is_zero(), (uRnd - uRnd) == 0
        (True, True)

        # Initialization from an element of Seqs.ore_algebra().base_ring()
        sage: Seqs(42), Seqs(42) == 42
        ([42, 42, 42, 42, ..., 42, ...], True)
        sage: u = Seqs(3*n**3 - n**2 + 5); u
        [5, 7, 25, 77, ..., 2111, ...]
        sage: (u[3], u[57], u[5:9])
        (77, 552335, [355, 617, 985, 1477])
        sage: try:
        ....:     u[5:2]
        ....: except IndexError as e:
        ....:     print(e)
        Upper index must not be smaller than the lower index
        sage: u[-5]

        # Extra conditions, getitem
        sage: fibo = Seqs([0,1], Sn**2 - Sn - 1); fibo
        [0, 1, 1, 2, 3, ..., 34, ...]
        sage: u = Seqs([1], Sn - 1); (u, u.is_one(), u == 1)
        ([1, 1, 1, ..., 1, ...], True, True)
        sage: v1 = Seqs([2], Sn-n); v2 = Seqs({0:2}, Sn-n); v3 = Seqs({1:2}, Sn-n)
        sage: v1 == v2, v1 == v3, v2 == v3
        (True, False, False)
        sage: v1,v2,v3
        ([2, 0, 0, 0, ..., 0, ...], [2, 0, 0, 0, ..., 0, ...], [2, 2, 4, ..., 725760, ...])
        sage: fibo2 = Seqs({1:1,2:1}, Sn**2-Sn-1); fibo2
        [1, 1, 2, 3, ..., 55, ...]
        sage: fibo == fibo2, fibo[1:9] == fibo2[1:9]
        (False, True)
        sage: fiboShift = Seqs([1,1], Sn**2-Sn-1); fiboShift
        [1, 1, 2, 3, ..., 55, ...]
        sage: fiboShift == fibo, fiboShift[0:9] == fibo[1:10]
        (False, True)
        sage: fiboAlt1 = Seqs({0:0,1:1,12:0}, Sn**2-Sn-1); fiboAlt1[12],fiboAlt1[0:12] == fibo[0:12]
        (0, True)
        sage: fiboAlt1[11:15]
        [89, 0, 89, 89]
        sage: fiboAlt2 = Seqs({0:0,1:1,12:0,18:0}, Sn**2-Sn-1); fiboAlt2[18], fiboAlt2[0:17] == fiboAlt1[0:17]
        (0, True)


        # Operations
        sage: fibo + 1
        [1, 2, 2, 3, ..., 35, ...]
        sage: eConsec = Seqs({1:1}, n*Sn-n-1); eConsec
        [1, 2, 3, ..., 10, ...]
        sage: fibo + eConsec
        [2, 3, 5, 7, ..., 65, ...]
        sage: eConsec - fibo
        [0, 1, 1, 1, ..., -45, ...]
        sage: eConsec * fibo # Note that the sequence starts at n=1 (because so does eConsec)
        [1, 2, 6, 12, ..., 306, ...]
        
    """
    Element = PRecursiveSequence
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
        if x in self.base_ring() :
            if x in self.base_ring().base_ring():
                # x is a constant
                return self.element_class (self, [x], self.generator()-1)
            # x is a polynomial
            gen = self._generator
            deg = x.degree()
            n = self.base_ring().gen()
            P = 0 
            # Proof of what follows : u(n) = P(n) => u(n+1) = P(n+1) <=> u(n+1) - u(n) - Q(n) = 0
            #       <=> u(n)*u(n+1) - u(n)*u(n) - Q(n)*u(n) = 0 <=> [P(n)*Sn - P(n) - Q(n)]*u(n) = 0
            for i in range(1,deg+1):
                for j in range(i):
                    P += x[i]*binomial(i,j)*n**j
            return self.element_class (self, [x.subs(0), x.subs(1)], x*gen - x - P)

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
        if self._ore_algebra.base_ring().has_coerce_map_from (S):
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
        # TODO : allow guessing with index that does not start at 1
        annihil = -guess (vals, self._ore_algebra)
        return self.element_class (self, vals, annihil)



