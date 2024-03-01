#! /home/porterll/python/Python-3.8.0/python

import sys
import numpy as np
import pickle
import shutil

def get_header(f):

    header = f[:2]

    block, i = get_block(f,2)
    header+=block

    return header, i

def get_block(f,idx):

    block = []

    while f[idx].split() != [] and idx < len(f)-1:
        block.append(f[idx])
        idx += 1

    return block, idx+1

if __name__ == '__main__':

    f = open(sys.argv[1]).read().splitlines()
    f2 = open(sys.argv[2],'rb')

    coevolved_alt = pickle.load(f2)

    header, sidx = get_header(f)

    blocks = []

    while sidx < len(f):
        block,sidx = get_block(f,sidx)
        blocks.append(block)

    ref_idxs = []

    for i in range(len(blocks)):
        info = blocks[i][0].split()
        ref_idxs += [(i,x) for x in range(len(info[1])) if info[1][x] != '-']

    try:
        toA = [ref_idxs[x] for x in coevolved_alt]
    except IndexError:
        shutil.copyfile(sys.argv[1],sys.argv[2].split('.')[0]+'/msas/'+sys.argv[1].split('/')[-1])
        sys.exit()
    
    relevant_blocks = np.unique(np.array([x[0] for x in toA]))

    batch = []
    for i in relevant_blocks:
        dat = []
        for j in range(1,len(blocks[i])):
            newstr = ''
            info = blocks[i][j].split()
            if info[0][:2] == '#=':
                continue
            dat.append([x for x in info[1]])

        dat = np.array(dat)
        for k in toA:
            if k[0] == i:
                try:
                    dat[:,k[1]] = 'A'
                except IndexError:
                    shutil.copyfile(sys.argv[1],sys.argv[2].split('.')[0]+'/msas/'+sys.argv[1].split('/')[-1])
                    sys.exit()

        batch.append(dat)

    for i in header:
        print(i)
    print('')

    space_found = 0
    indent = None
    rb_idx = 0

    for i in range(len(blocks[0][1])):
        if space_found and blocks[0][1][i] != ' ':
            indent = i
            break
        elif blocks[0][1][i] == ' ':
            space_found = 1


    for i in range(len(blocks)):
        if i not in relevant_blocks:
            for j in blocks[i]:
                print(j)
        else:
            sidx = 0
            print(blocks[i][0])
            for j in range(1,len(blocks[i])):
                if blocks[i][j][0] == '#':
                    print(blocks[i][j])
                else:
                    print(blocks[i][j][:indent]+''.join(batch[rb_idx][sidx]))
                    sidx += 1
            rb_idx+=1

        if i != len(blocks)-1:
            print('')
        else:
            print('//')
        





        
