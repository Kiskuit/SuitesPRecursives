from __future__ import print_function

from seq import *
from parentSeq import *
from ore_algebra import *
import time


def test1():
    R,n = ZZ['n'].objgen()
    A,Sn = OreAlgebra(R, 'Sn').objgen()

    a1 = Sn**2 - Sn - 1
    init1 = [0,1]

    range_ = 1000

    sum_ = 0
    for _ in range (range_):
        begin = time.time()
        m1 = a1.to_list (init1, 100)[-1]
        end = time.time()
        sum_ += end - begin
    print ("Single call to to_list() :",sum_)
    
    sum_ = 0
    for _ in range(range_):
        begin = time.time()
        init1 = [0,1]
        for i in range(2,100):
            a1.to_list (init1, 3, start=i-2,append=True)
            init1 = init1[1:]
        m2 = init1[-1]
        end = time.time()
        sum_ += end - begin
    print ("Multiple calls to to_list() :", sum_)
    if m1!=m2:
        print("ERROR")
        print(m1, "\n", m2)
    else :
        print ("SUCCESS")

def test2():
    pass

if __name__ == "__main__":
    #n = ZZ['n'].gen()
    #P = ParentSeqRec (ZZ['n'], 'Sn', RR)
    #a = P ([0,1], Sn^2 - Sn - 1)

    test2()
    
