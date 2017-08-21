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
###   import urllib2
###   from bs4 import BeautifulSoup
import re
import csv
import parent_p_recursive_sequences as ppseq
# Sage imports
from sage.rings.integer_ring import ZZ
import ore_algebra

# I hope the html is build automatically and the tables[4] is always the one of interest
# TODO use tables[3] instead coz we need the index
magicNumber = 4
maxSeq = 289480

def getSequenceFromTxt(seq_filename, offset_filename):
    dbfile = "test.txt"
    seq_fd = open (seq_filename, 'r')
    offset_fd = open(offset_filename, 'r')
    # TODO iterate though the min of both files?
    for line in seq_fd:
        seq = parseSeq(seq_fd.readline())
        start = parseOffset(offset_rd.readline())
        annihilator = guess(seq, start)
        if annihilator is not None :
            order = annihilator.order()
            cond = seq[:order]
            with open(dbfile, 'a') as fd:
                fd.write("A%06d,%s,%s\n" % (indexSeq, annihilator, cond))

def parseOffset(line):
    line = line[11:]
    line = line.split(",")[0]
    return int(line)

def parseSeq(line):
    # TODO
    return [0,1,1,2,3,5,8,13,21,34,55]

def getSequenceFromPage(url):
    # Open url and get content
    sloanePage = urllib2.urlopen(url)
    sloaneSoup = BeautifulSoup(sloanePage, "lxml")
    # Find the right table
    tables = sloaneSoup.find_all("table")
    # table = tables[magicNumber]
    table = tables[3]
    table = table.find_all("tt")[2:]
    # Remove whitespace and brackets
    pattern = re.compile(r'[^0-9,]')
    table = map(lambda x: int(re.sub(pattern,'',x.string)), table)
    seq = {}
    k,v = table[0::2], table[1::2]
    if k[-1] != k[0] + len(k) -1:
        # Case where the indices are not consecutives
        return None
    return (k[0],v)

def browsePages(start, stop):
    dbfile = '../db/properties.csv'
    # Progress vars
    i_progress = 1
    progress = int(stop-start)/int(10)
    # Empty file
    open(dbfile,'w').close()
    # Browse all pages
    for indexSeq in range(start, stop):
        url = "https://oeis.org/A%06d/list" % indexSeq
        if indexSeq == i_progress*progress:
            print ("progress : %d0%%" % i_progress)
            i_progress += 1
        seq = getSequenceFromPage(url)
        if seq is None:
            continue
        else:
            start,seq = seq
        annihilator = guess(seq, start)
        if annihilator is not None :
            order = annihilator.order()
            cond = seq[:order]
            with open(dbfile, 'a') as fd:
                fd.write("A%06d,%s,%s\n" % (indexSeq, annihilator, cond))


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

###   if __name__ == "__main__":
###       #print (getSequenceFromPage("https://oeis.org/A005206/list"))
###       from sys import argv
###       if len(argv) == 1:
###           start = 1
###           stop = maxSeq
###       elif len(argv) == 2:
###           start = 1
###           # TODO
###           stop = 200
###           ###   stop = int(argv[1])
###       ###   elif len(argv) > 2:
###       ###       start = int(argv[1])
###       ###       stop = int(argv[2])
###       else:
###           ValueError("Please call this script ./database.py `start` `stop` or ./database `stop`")
###       if stop < start or start < 0:
###           raise ValueError("Please provide positive integers such that start < stop")
###       browsePages(start, stop)
    
if __name__ == "__main__":
    getSequenceFromTxt("seq.txt", "offset.txt")

