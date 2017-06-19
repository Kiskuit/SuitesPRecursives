from __future__ import print_function

import time
import csv
from sage.all import *
from ore_algebra import *

if __name__ == "__main__":
    # TODO What is 2nd param of forward_matrix for?
    A,n = ZZ["n"].objgen()
    R,Sn = OreAlgebra(A,"Sn").objgen()

    fibo = Sn**2 - Sn - 1

    maxOrder = 5
    minOrder = 1
    fiboTime = []
    execTime = [[] for _ in range(minOrder,maxOrder+1)]
    cond = [range(i) for i in range(minOrder,maxOrder+1)]
    for i in [10, 100, 1000]:
        lst = 0
        fwdMtx = 0
        nbExec = 10
        for k in range(nbExec):
            startLst = time.time()
            fibo.to_list([0,1],i)
            stopLst = time.time()
            
            startMtx = time.time()
            P,Q = fibo.forward_matrix_bsplit (i,0)
            P*Matrix([[0],[1]])/Q
            stopMtx = time.time()

            lst += (stopLst - startLst)*1000./nbExec
            fwdMtx += (stopMtx - startMtx)*1000./nbExec

        fiboTime.append((lst, fwdMtx))
        print (i, "Done :", lst, "\t", fwdMtx)

        for j in range(minOrder,maxOrder+1):
            k = 0
            while k<nbExec :
                f = R.random_element (j)
                # Using to_list method
                startLst = time.time()
                _ = f.to_list (cond[j-minOrder], i)
                stopLst = time.time()

                # Using forward_matrix
                startMtx = time.time()
                P,Q = f.forward_matrix_bsplit (i,0)
                if Q==0:
                   continue
                _ = P*Matrix([[e] for e in cond[j-minOrder]])/Q
                stopMtx = time.time()

                lst += (stopLst - startLst)*1000./nbExec
                fwdMtx += (stopMtx - startMtx)*1000./nbExec
                k += 1
            execTime[j-minOrder].append((lst, fwdMtx))
            print (i, "Done :", lst, "\t", fwdMtx)


    with open('bench.csv','w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        times = zip(*execTime)
        for i,j in enumerate([10,100,1000]):
            lstTime, mtxTime = zip(*times[i])
            writer.writerow ([j] + ['to_list'] + [fiboTime[i][0]] + list(lstTime))
            writer.writerow ([j] + ['forward_matrix'] + [fiboTime[i][1]] + list(mtxTime))
