#! /usr/local/bin/python3
#
#
#This faster implementation of SPEACH-AF (Stein and Mchaourab, PLOS Comp. Bio. 2022)
#finds all unique sets of alanines to knock out in an input PDB file.
#It works by moving down a protein chain in steps of one, finding all contacts within an 11-residue
#window, W (excluding contacts +/- 4 residues adjacent to W).
#Contacts are defined as atoms with 4 Angstroms of any atom in W.
#Unique sets of contacts are saved as pickled files within a user-specified directory.

#Author: Lauren Porter (porterll@nih.gov)

import sys
import numpy as np
from scipy.spatial import KDTree
import pandas as pd
import pickle
import os

#Molecule module taken from Srinivasan and Rose, PNAS 1999
class Molecule:
    """Class to describe a molecule.

    An instance of Molecule has the following attibutes:

        o *filename* - name of the file from which the molecule was read

        o *name* - name/title for the molecule (string) (default '')

        o *numatm* - number of atoms in the molecule (integer)

        o *numres* - number of residues in the molecule (integer)

        o *resptr* - numpy array of type 'i' and length *numres*. This
        array contains the index of the first atom of each residue in the
        molecule

        o *atmidx* - numpy array of type 'i' and length *numatm*. This
        array contains the index of each atom in the molecule as specified
        in the PDB file

        o *atmnam* - A python list of length *numatm*.  Contains the name
        of each atom in the molecule as a string of length 4, following the
        PDB convention

        o *resnam* - A python list of length *numres*.  Contains the name
        of each residue in the molecule

        o *resnum* - A numpy array of type 'i' and length *numatm*. Contains
        the residue number of each atom in the molecule as specified in the
        PDB file

        o *coords* - A numpy matrix of type 'f' and shape (*numatm*, 3).
        Contains the cartesian coordinates of each atom in the molecule.

    """

    def __init__(self, filename):
        """Create a new Molecule instance.

        Arguments

            o filename - name of file containing the coordinates of the
                         protein in PDB format

        """

        if type(filename) == type(''):
            self.filename = os.path.abspath(os.path.expanduser(
                os.path.normpath(filename)))

            if self.filename[-3:] == '.gz':
                import gzip
                myOpen = gzip.open
            else:
                myOpen = open

            f = myOpen(self.filename)
            data = f.readlines()
            f.close()

        else:
            data = filename.readlines()
            self.filename = 'data.pdb'

        # count up number of atoms and residues in molecule.
        # also store the molecule name and the indices of the
        # first atom of each residue

        self.name = ''
        self.trial = 0

        self.numatm = self.numres = 0

        self.resptr = []
        self.pdbidx = []
        self.pdbicd = []
        self.pdbchn = []
        
        oldres = None
        oldchn = None

        for line in data:
            if line[:6] == 'COMPND':
                self.name = str.strip(line[6:])
            elif line[:4] == 'ATOM':
                if oldchn != None and line[21] != oldchn:
                    break
                if line[22:27] != oldres:
                    self.numres = self.numres + 1
                    self.resptr.append(self.numatm)
                    self.pdbidx.append(int(line[22:26]))
                    self.pdbicd.append(line[26])
                    self.pdbchn.append(line[21])
                    oldres = line[22:27]
                    oldchn = line[21]
                self.numatm = self.numatm + 1

        self.resptr.append(self.numatm)
        self.pdbidx = np.array(self.pdbidx, 'i')
        self.resptr = np.array(self.resptr, 'i')
        
        # allocate all fields for Molecule structure

        # -- atom indices as specified by the pdb file

        self.atmidx = np.zeros(self.numatm, 'i')

        # -- atom names

        self.atmnam = ['']*self.numatm

        # -- residue names

        self.resnam = ['']*self.numres

        # -- residue numbers

        self.resnum = np.zeros(self.numatm, 'i')

        # -- x, y, z coordinates

        self.coords = np.zeros((self.numatm, 3), 'f')

        # parse the data

        self._parse(data)

    def _parse(self, data):

        na = 0
        nr = -1
        oldres = None
        oldchn = None

        for line in data:
            if line[:4] == 'ATOM':
                if oldchn != None and line[21] != oldchn:
                    if verbose:
                        print(('Using chain %s, ignoring all others' % oldchn))
                    break
                oldchn = line[21]
                
                rnam = str.strip(line[17:21])
                anam = self.atmnam[na] = line[12:16]
                self.coords[na] = [float(line[30:38]),
                                   float(line[38:46]),
                                   float(line[46:54])]
                
                if line[22:27] != oldres:
                    nr = nr + 1
                    self.resnam[nr] = rnam
                    oldres = line[22:27]

                self.resnum[na] = nr
                    
                na = na + 1

            elif line[:6] == 'COMPND':
                self.name = line[6:].strip()


#Find all contacts in an 11-residue window.  If unique, save as a picke file in the user-specified directory
def get_contact_lists(name,cres,indices,m,contact_sets):

    cidx = cres
    idxs2avoid = range(m.resptr[max(0,cidx-9)],m.resptr[min(m.numres,cidx+10)])
    contacts = []
 
    for i in range(m.resptr[max(0,cidx-5)],m.resptr[min(cidx+6,m.numres)]):
        idxs = [m.resnum[x] for x in indices[i] if x not in idxs2avoid]
        new = np.unique(np.array([x for x in idxs if x not in contacts]))
        contacts += list(new)
 
    contacts.sort()

    if contacts and contacts not in contact_sets:
        contact_sets.append(contacts)
        print(sys.argv[3]+'/'+name+'_%i_%i.pik'%(max(0,cidx-5),min(cidx+6,m.numres)))
        print(contacts)
        of = open(sys.argv[3]+'/'+name+'/'+name+'_%i_%i.pik'%(max(0,cidx-5),min(cidx+6,m.numres)),'wb')
        pickle.dump(contacts,of)
        of.close()
    
    return contact_sets

#Find all contacts within a PDB file
def get_contacts(name,minidx,maxidx):
    
    m = Molecule(name+'_clean.pdb')
    kdXYZ = KDTree(m.coords)
    indices = kdXYZ.query_ball_tree(kdXYZ, 4.0) #Any pair of atoms <= 4.0 Angstroms apart are considered in contact

    contact_sets = []

    os.mkdir(sys.argv[3]+'/'+name.split('/')[-1][:4])
    
    for i in range(minidx,maxidx):
        contact_sets = get_contact_lists(name.split('/')[-1][:4],i,indices,m,contact_sets)
 
if __name__ == '__main__':

    if len(sys.argv) < 4:
        sys.exit('Usage: ./get_contact_sets_noLINUS.py <input file with list of pdbs> <path to pdbs> <output path>')
 
    fs_pairs = open(sys.argv[1]).read().splitlines()
    prefix = sys.argv[2]

    for i in range(1,len(fs_pairs)):
        info = fs_pairs[i].split(',')

        get_contacts(prefix+'/'+info[0][:4],int(info[-2].split('-')[0])-1,int(info[-2].split('-')[1])-1)
        get_contacts(prefix+'/'+info[1][:4],int(info[-1].split('-')[0])-1,int(info[-1].split('-')[1])-1)



