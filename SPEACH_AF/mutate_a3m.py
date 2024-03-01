#! /home/porterll/python/Python-3.8.0/python

import sys
import numpy as np
import pickle

if __name__ == '__main__':

    f = open(sys.argv[1]).read().splitlines()

    f2 = open(sys.argv[2],'rb')
    coevolved_alt = pickle.load(f2)

    print(f[0])
    print(f[1])

    for i in range(2,len(f)):
        idx = 0
        fstr = ''
        if f[i][0] == '>':
            print(f[i])
        else:
            for j in range(len(f[i])):
                if f[i][j] == '-' or f[i][j].isupper():
                    if idx in coevolved_alt:
                        fstr += 'A'
                    else:
                        fstr += f[i][j]
                    idx += 1

                else:
                    fstr += f[i][j]
            print(fstr)
                
