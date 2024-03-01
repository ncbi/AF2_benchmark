import sys
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("name", help="name to save everything under")
parser.add_argument("--target_list", nargs='*', help="List of target names to run")
parser.add_argument("--targets_file", default="", help="File with list of target names to run")
parser.add_argument("--recycles", type=int, default=1, help="Number of recycles when predicting")
parser.add_argument("--model_num", type=int, default=1, help="Which AF2 model to use")
parser.add_argument("--seed", type=int, default=0, help="RNG Seed")
parser.add_argument("--verbose", action='store_true', help="print extra")
parser.add_argument("--deterministic", action='store_true', help="make all data processing deterministic (no masking, etc.)")
parser.add_argument("--use_native", action='store_true', help="add the native structure as a decoy, and compare outputs against it")
parser.add_argument("--mask_sidechains", action='store_true', help="mask out sidechain atoms except for C-Beta")
parser.add_argument("--mask_sidechains_add_cb", action='store_true', help="mask out sidechain atoms except for C-Beta, and add C-Beta to glycines")
parser.add_argument("--seq_replacement", default='', help="Amino acid residue to fill the decoy sequence with. Default keeps target sequence")
parser.add_argument("--af2_dir", default="/n/home01/jroney/alphafold-latest/", help="AlphaFold code and weights directory")
parser.add_argument("--decoy_dir", default="/n/home01/jroney/template_injection/decoy_set/", help="Rosetta decoy directory")
parser.add_argument("--output_dir", default="/n/ovchinnikov_lab/Lab/af2rank/outputs/", help="Rosetta decoy directory")
parser.add_argument("--tm_exec", default="/n/home01/jroney/tmscore/TMscore", help="TMScore executable")

args = parser.parse_args()

sys.path.append(args.af2_dir)

from alphafold.model import model
from alphafold.model import config
from alphafold.model import data

from alphafold.data import parsers
from alphafold.data import pipeline

from alphafold.common import protein
from alphafold.common import residue_constants


import jax
import jax.numpy as jnp
import numpy as np
import re
import subprocess
from collections import namedtuple

'''import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--af2_dir", default="/n/home01/jroney/alphafold-latest/", help="AlphaFold code and weights directory")

args = parser.parse_args()

sys.path.append(args.af2_dir)

from alphafold.model import model'''


