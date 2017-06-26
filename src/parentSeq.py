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
from sage.all import *
from ore_algebra import OreAlgebra, guess
from seq import PRecSequence

# TODO change class it inherits from ?
class ParentSeqRec (Ring, UniqueRepresentation):
    Element = PRecSequence
    # TODO generator as optional arg and if none : gen = 'S'+base.gen() ?
    def __init__(self, base_ring, generator, val_ring, category=None):
        self._val_ring = val_ring
        # TODO catch error, to raise mine
        self._ore_algebra, self._generator = OreAlgebra (base_ring, generator).objgen()
        Ring.__init__(self, base_ring, category=category or Rings())

    def val_ring(self):
        return self._val_ring

    def ore_algebra(self):
        return self._ore_algebra

    def generator(self):
        return self._generator

    def base_ring(self):
        return self.ore_algebra().base_ring()

    def _repr_(self):
        ret = "P-Recursive sequences over `{}`, with values in `{}`".format(self.base_ring(),
                                                                            self.val_ring())
        return ret

    def _element_constructor_(self, *args, **kwargs):
        if len(args) != 1:
            return self.element_class (self, *args, **kwargs)
        x = args[0]
        try :
            P = x.parent()
        except AttributeError :
            # TODO what if attributerror?
            return self.element_class (self, x, **kwargs)

        # x is in base_ring()
        # TODO take into account cases where there are roots to the pol
        # TODO lower boundary guessing --> `8` hard coded but ugly
        if x in self.base_ring() :
            gen = self.base_ring().gen()
            condInit = x.subs({gen:0})
            annihil = -guess ([x.subs({gen:i}) for i in range(8)], self._ore_algebra)
            return self.element_class (self, condInit, annihil)

        # x is a constant of val_ring
        #   return with annihilator Sn - 1 (where Sn is replaced by actual generator)
        if x in self._val_ring :
            return self.element_class (self, {0:x}, self._ore_algebra.gen() - 1)

        # Default case
        return self.element_class (self, x, **kwargs)

    def _coerce_map_from_ (self, S):
        if self.ore_algebra().has_coerce_map_from (S):
            return True
        if self._val_ring.has_coerce_map_from (S):
            return True
        
        # Default case 
        return False



if __name__ == "__main__":
    R = ZZ['n']
    par = ParentSeqRec (R,'Sn',RR)
    print ("-"*32)
    print (par)
    print ("-"*32)
    try:
        par = ParentSeqRec (R,'Sx',RR)
        print (par)
    except TypeError:
        print ("Sx is not a generator of the OreAlgebra over {}.".format(R))
    print ("-"*32)
    try:
        par = ParentSeqRec (ZZ, 'Sn', RR)
        print (par)
    except TypeError:
        print ("ZZ isnt a proper ring for OreAlgebra.")
    print ("-"*32)
