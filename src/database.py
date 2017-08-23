#!/usr/bin/python2
# -*- coding: utf8 -*-

# ---------------------------------------------------------
# Author : Mathis Caristan (mathis.caristan@gmail.com)
# Date : 12/07/2017
# Class : database
# Copyright 2017 Mathis Caristan
# ---------------------------------------------------------
from __future__ import print_function

# Python imports
import parent_p_recursive_sequences as ppseq
# Sage imports
from sage.rings.integer_ring import ZZ
import ore_algebra

def parseFiles(seq_filename, offset_filename, dbfile="out.txt", maxSeq=999999):
    open(dbfile,"w").close()
    seq_fd = open (seq_filename, 'r')
    offset_fd = open(offset_filename, 'r')
    # Iterates through lines of the offset file, as it should be the shortest
    for indexSeq,line in enumerate(offset_fd):
        if indexSeq > maxSeq:
            break
        seq = parseSeq(seq_fd.readline())
        start = parseOffset(line)
        annihilator = guess(seq, start)
        if annihilator is not None :
            order = annihilator.order()
            cond = seq[:order]
            with open(dbfile, 'a') as fd:
                fd.write("A%06d,%s,%s\n" % (indexSeq+1, annihilator, cond))
    seq_fd.close()
    offset_fd.close()

def parseOffset(line):
    ret = line.split(" ")[-1]
    ret = ret.split(",")[0]
    return int(ret)

def parseSeq(line):
    ret = filter(lambda x: x is not "", line.strip().split(','))[1:]
    return map(lambda x: int(x), ret)

def guess(seq, start):
    Seqs = ppseq.ParentPRecursiveSequences(ZZ['n'])
    
    for i in range(10,30):
        alg = Seqs.ore_algebra()
        if i > len(seq):
            return None # TODO handle it better?
        try:
            annihilator = ore_algebra.guess(seq[:i], alg, start=start)
        except ValueError:
            continue
        # TODO check if annihilator is good or not
        seq_guessed = annihilator.to_list(seq[:annihilator.order()], len(seq), start=start)
        if seq_guessed == seq :
            return annihilator
    return None

if __name__ == "__main__":
    # from sys import argv
    # TODO parse arguments
    parseFiles("../db/seq.txt", "../db/offset.txt", maxSeq=200)


