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
import urllib2
from bs4 import BeautifulSoup
import re
import csv
import parent_p_recursive_sequences as ppseq
# Sage imports
from sage.rings.integer_ring import ZZ
import ore_algebra

# For tests
import random

# I hope the html is build automatically and the tables[4] is always the one of interest
# TODO use tables[3] instead coz we need the index
magicNumber = 4
maxSeq = 201006

def getSequenceFromPage(url):

    # Open url and get content
    sloanePage = urllib2.urlopen(url)
    sloaneSoup = BeautifulSoup(sloanePage, "lxml")
    # Find the right table
    tables = sloaneSoup.find_all("table")
    table = tables[magicNumber]
    # Remove whitespace and brackets
    pattern = re.compile(r'[^0-9,]')
    string = re.sub(pattern, '', table.pre.string)
    # Turn the string into an array
    seq = string.split(",")
    seq = map(lambda x:int(x),seq)

    return seq

def browsePages():
    maxSeq = 200
    dbfile = '../db/properties.csv'
    # Progress vars
    i_progress = 1
    progress = int(maxSeq)/int(10)
    # Empty file
    open(dbfile,'w').close()
    # Browse all pages
    for indexSeq in range(1,maxSeq):
        url = "https://oeis.org/A%06d/list" % indexSeq
        if indexSeq == i_progress*progress:
            print ("progress : %d0%%" % i_progress)
            i_progress += 1
        seq = getSequenceFromPage(url)
        
        annihilator = guess(seq)
        if annihilator is not None :
            order = annihilator.order()
            cond = seq[:order]
            with open(dbfile, 'a') as fd:
                fd.write("A%06d,%s,%s\n" % (indexSeq, annihilator, cond))


def guess(seq):
    Seqs = ppseq.ParentPRecursiveSequences(ZZ['n'])
    
    # if random.randint(1,5) == 1: # Simulation of guessing successful or not
    #     return "Sn**2-Sn-1"
    for i in range(10,30):
        alg = Seqs.ore_algebra()
        if i > len(seq):
            return None # TODO handle it better?
        try:
            annihilator = ore_algebra.guess(seq[:i], alg)
        except ValueError:
            continue
        # TODO check if annihilator is good or not
        seq_guessed = annihilator.to_list(seq[:annihilator.order()], len(seq))
        if seq_guessed == seq :
            return annihilator
    return None

if __name__ == "__main__":
    # getSequenceFromPage("https://oeis.org/A005206/list")
    browsePages()
    
