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
from ore_algebra import OreAlgebra
from seq import PRecSequence

# TODO change class it inherits from ?
class ParentSeqRec (Ring):
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
