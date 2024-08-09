import glob
import itertools
import numpy as np

import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, Selection
from scipy.spatial import distance_matrix

class CMAP():
    """
    Take in a PDB structure that was loaded using BioPython and produce a contact map
    """

    def __init__(self,structure):
        target_file = structure; pdb = structure.split('/')[-1].replace('.pdb', '')
        structure = PDBParser().get_structure(pdb, target_file)
        self.structure = structure


    def Create_Map(self, cutoff):
        """
        Create the contact map from the PDB structure
        gaps should be filled in and homo-oligomer chain lengths should match
        even if crystal density was not present
        """
        chains = {}
        for chain in self.structure[0].get_chains():
            residues = Selection.unfold_entities(chain, 'R')
            residues = sorted(residues, key=lambda r:r.get_id()[1])
            chains[chain] = [residue['CA'] for residue in residues if 'CA' in residue.child_dict.keys()]

        coor_sep = np.array([[atom.coord for atom in chain[1]] for chain in chains.items()],dtype=object)
        coor_idx = [[i for i,a in enumerate(chain)] for chain in coor_sep]
        tot_idx,offset = [],0
        for idx, coor in enumerate(coor_idx):
            if idx == 0:
                tot_idx.append(list(coor_idx[idx]))
            else:
                offset = offset + coor_idx[idx - 1][-1] +1
                tot_idx.append([pos+offset for pos in coor_idx[idx]])
        tot_idx_flat = [pos for chain in tot_idx for pos in chain]

        coor = np.array([pos for chain in coor_sep for pos in chain])
        dist_matrix = distance_matrix(coor[tot_idx_flat,:],coor[tot_idx_flat,:],p=2)

        #contact map
        contact_map = np.zeros((len(tot_idx_flat),len(tot_idx_flat)))
        contact_map[dist_matrix < cutoff] = 1
        contact_map[dist_matrix == 0] = 0

        res_list = [res for key in chains for res in chains[key]]
        combo = list(itertools.permutations(res_list,2))
        res_idx = list(itertools.permutations(list(range(len(res_list))),2))

        fill = []
        for key,idx in zip(combo,res_idx):
            key_a,key_b = key[0],key[1]
            if key_a.parent is not None and key_b.parent is not None and idx[1] > idx[0]:
                a = np.array([key_a.parent.child_dict[key].coord for key in key_a.parent.child_dict])
                b = np.array([key_b.parent.child_dict[key].coord for key in key_b.parent.child_dict])
                dist_temp = distance_matrix(a,b,p=2)
                if np.sum(dist_temp <= cutoff).astype(bool):
                    fill.append(idx)
        for i,j in fill:
            contact_map[i][j] = 1

        #last remove hits within +-3 of the diagonal
        mask = np.zeros((len(tot_idx_flat),len(tot_idx_flat)))
        mask = np.abs(np.arange(len(tot_idx_flat)) - np.arange(len(tot_idx_flat))[:,np.newaxis]) <= 3
        contact_map[mask] = 0

        #upper triangle is a contact map that represents every contact between any heavy atom of two residues
        #lower triangle is a contact map that represents contacts between CAs of two residues
        tu = np.triu_indices(contact_map.shape[0])
        contact_map[tu[::-1]] = contact_map[tu]

        #create a mask that makes oligomers easier to visualize
        #as is this will only handle up to 5 chains, makes interchain squares darker
        mask = np.zeros((len(tot_idx_flat),len(tot_idx_flat)))
        perm = list(range(len(tot_idx)))
        perm = [p for p in itertools.permutations(perm, r=2) if p[0] < p[1]]

        value = 0.2
        for i in perm:
            mask_idx = np.array(list(itertools.product(tot_idx[i[0]],tot_idx[i[1]])))

            mask[mask_idx[:,0],mask_idx[:,1]] = 1 - int(i[1]-i[0])*value
            mask[mask_idx[:,1],mask_idx[:,0]] = 1 - int(i[1]-i[0])*value

        return contact_map,mask

def pad_with(vector,pad_width,iaxis,kwargs):
    """
    custom paramter to create a padded edge of a given value around an np.array using np.pad()
    """
    pad_value = kwargs.get('padder', 0)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value

def overlap_definition(A, B, mtx_return=False):
    #pad edges for sliding window consistency
    padded_A = np.pad(A, ((1,1),(1,1)), mode='constant', constant_values=0)
    #Extract all possible windows
    windows = np.lib.stride_tricks.sliding_window_view(padded_A, (3,3))
    #positions with contacts
    mask_1 = np.where(B == 1)
    mask_2 = np.where(B == 2)
    mask_3 = np.where(B == 3)
    #find number of contacts in B that match a (3,3) window in A that has a nonzero elements
    matches_1 = np.any(windows[mask_1] != 0, axis=(-1,-2))
    matches_2 = np.any(windows[mask_2] != 0, axis=(-1,-2))
    matches_3 = np.any(windows[mask_3] != 0, axis=(-1,-2))

    # return B array with +s where a +-1 match in A exists and -s where they don't exist
    # 0s where no contact exists, this function retains the type of contact
    if mtx_return == True:
        B_true = np.copy(B)
        B_true[mask_1] = matches_1 * 1
        B_true[mask_2] = matches_2 * 2
        B_true[mask_3] = matches_3 * 3

        B_false = np.copy(B)
        B_false[mask_1] = ~matches_1 * 1
        B_false[mask_2] = ~matches_2 * 2
        B_false[mask_3] = ~matches_3 * 3
        B_false = B_false*-1
        return B_true + B_false
    return sum([*matches_1, *matches_2, *matches_3])


def MTX_comparison(cmap_down, cmap_up, mtx):
    up_rgba = np.full((*cmap_up.shape, 4), (255/255, 168/255, 0/255, 0) , dtype='f') #orange
    up_rgba[np.where(cmap_up == 1)] = (255/255, 168/255, 0/255, 1)

    down_rgba = np.full((*cmap_down.shape, 4), (0/255, 87/255, 255/255, 0) , dtype='f') #blue2
    down_rgba[np.where(cmap_down == 1)] = (0/255, 87/255, 255/255, 1)

    tri_upper_idx = np.triu_indices(mtx.shape[0])
    tri_lower_idx = np.tril_indices(mtx.shape[0])

    up_rgba[np.where((up_rgba[:,:,-1] == 1) & (down_rgba[:,:,-1] == 1))] = (0,0,0,1)
    down_rgba[np.where((up_rgba[:,:,-1] == 1) & (down_rgba[:,:,-1] == 1))] = (0,0,0,1)

    rgba_cmap = np.empty_like(up_rgba)
    rgba_cmap[tri_lower_idx] = up_rgba[tri_lower_idx]
    rgba_cmap[tri_upper_idx] = down_rgba[tri_upper_idx]

    up_rgba[tri_lower_idx[0], tri_lower_idx[1], -1] = up_rgba[tri_lower_idx[0], tri_lower_idx[1], -1] * mtx[tri_lower_idx]
    down_rgba[tri_upper_idx[0], tri_upper_idx[1], -1] = down_rgba[tri_upper_idx[0], tri_upper_idx[1], -1] * mtx[tri_upper_idx]

    rgba_zscore = np.empty_like(up_rgba)
    rgba_zscore[tri_lower_idx] = up_rgba[tri_lower_idx]
    rgba_zscore[tri_upper_idx] = down_rgba[tri_upper_idx]

    return rgba_cmap, rgba_zscore

def OCTO_PLOT(cmap_down, cmap_up, msa_tr, distogram, name):

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(18,8))

    for idx, (msa, dist) in enumerate(zip(msa_tr, distogram)):
        msa = np.load(msa)
        dist = np.load(dist)
        rgba_cmap, rgba_zscore = MTX_comparison(cmap_down, cmap_up, msa) 

        #top
        ax[0][idx].imshow(rgba_zscore, interpolation='none')
        ax[0][idx].invert_yaxis()

        #bottom
        mask = np.abs(np.arange(dist.shape[0]) - np.arange(dist.shape[0])[:,None]) <= 2
        dist[mask] = np.nan

        ax[1][idx].imshow(dist, interpolation='none', cmap = 'Greys')
        ax[1][idx].imshow(rgba_cmap, interpolation='none', alpha = 0.35)
        ax[1][idx].invert_yaxis()
    plt.savefig(f'{name}.png', dpi = 600)
    plt.clf()
    plt.close()



dir_path = glob.glob("rfah/*")
numpy_files = [file for file in dir_path if '.npy' in file]

dist_files = sorted([file for file in dir_path if 'distogram' in file and '.npy' in file])
coev_files = sorted([file for file in dir_path if 'extra' in file and '.npy' in file])

#rfah structures
rfah_beta  = 'rfah/RfaH_unrelaxed_rank_005_alphafold2_model_1_seed_001.r3.pdb'
rfah_alpha = 'rfah/RfaH_unrelaxed_rank_004_alphafold2_model_2_seed_000.r3.pdb'
down_structure = rfah_beta; up_structure = rfah_alpha;
states = ['active', 'autoinhibited']

Cmap_up = CMAP(up_structure) #cmap object now contains structural information!!!
cmap_up,mask_up = Cmap_up.Create_Map(cutoff=8) # cutoff in angstom
Cmap_down = CMAP(down_structure) #cmap object now contains structural information!!!
cmap_down,mask_down = Cmap_down.Create_Map(cutoff=8) # cutoff in angstom

OCTO_PLOT(cmap_down, cmap_up, msa_tr=coev_files, distogram=dist_files, name='rfah')


# no recycles
dir_path = glob.glob("rfah_no_recycles/*")
numpy_files = [file for file in dir_path if '.npy' in file]

dist_files = sorted([file for file in dir_path if 'distogram' in file and '.npy' in file])
coev_files = sorted([file for file in dir_path if 'extra' in file and '.npy' in file])

OCTO_PLOT(cmap_down, cmap_up, msa_tr=coev_files, distogram=dist_files, name='rfah_no_recycles')
