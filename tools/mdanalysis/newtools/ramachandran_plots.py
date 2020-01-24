#!/usr/bin/env python

import argparse
import csv
import sys
import yaml
from collections import namedtuple

import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals

import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt

import numpy as np

import seaborn as sns


def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--itraj', help='input traj')
    parser.add_argument('--istr', help='input str')
    parser.add_argument('--itrajext', help='input traj ext')
    parser.add_argument('--istrext', help='input str ext')
#    parser.add_argument('--isegid1', help='segid 1')
#    parser.add_argument('--iresid1', help='resid 1')
#    parser.add_argument('--iname1', help='name 1')
#    parser.add_argument('--isegid2', help='segid 2')
#    parser.add_argument('--iresid2', help='resid 2')
#    parser.add_argument('--iname2', help='name 2')
#    parser.add_argument('--isegid3', help='segid 3')
#    parser.add_argument('--iresid3', help='resid 3')
#    parser.add_argument('--iname3', help='name 3')
#    parser.add_argument('--isegid4', help='segid 4')
#    parser.add_argument('--iresid4', help='resid 4')
#    parser.add_argument('--iname4', help='name 4')
#    parser.add_argument('--isegid5', help='segid 1')
#    parser.add_argument('--iresid5', help='resid 1')
#    parser.add_argument('--iname5', help='name 1')
#    parser.add_argument('--isegid6', help='segid 2')
#    parser.add_argument('--iresid6', help='resid 2')
#    parser.add_argument('--iname6', help='name 2')
#    parser.add_argument('--isegid7', help='segid 3')
#    parser.add_argument('--iresid7', help='resid 3')
#    parser.add_argument('--iname7', help='name 3')
#    parser.add_argument('--isegid8', help='segid 4')
#    parser.add_argument('--iresid8', help='resid 4')
#    parser.add_argument('--iname8', help='name 4')
    parser.add_argument('--iyml', help='input in yml format')
    parser.add_argument('--output', help='output')
    parser.add_argument('--oramachandran_plot', help='dihedral plot')
    return parser.parse_args()


args = parse_command_line(sys.argv)

import yaml
with open(args.iyml) as file:
    params = yaml.load(file, Loader=yaml.FullLoader)

Dihedral = namedtuple(
    'Dihedral', ['atom1', 'atom2', 'atom3', 'atom4'])

for k,v in params.items():
   for a in ['phi','psi']:
       assert (a in v), "Key %s is missing in inputs: %s " % (a, k)
       atoms=[]
       for b in ['atom1','atom2','atom3','atom4']:
          assert (b in v[a]), "Key %s is missing in inputs: %s %s" % (b, k, a)
          for c in ['segid','resid','name']:
              assert (c in v[a][b]), "Key %s is missing in inputs: %s %s %s " % (c, k, a, b)
          atoms.append("(segid %s and resid %s and name %s)" % (v[a][b]['segid'], v[a][b]['resid'], v[a][b]['name']))
       print(atoms)
       if a=='phi':
           dihe_phi = Dihedral(atoms[0],atoms[1], atoms[2], atoms[3])
       if a=='psi':
           dihe_psi = Dihedral(atoms[0],atoms[1], atoms[2], atoms[3])

# order of dihedral atom is the crystallographic definition
# (see glycanstructure.org)

assert(dihe_phi), "phi dihedral doesn't exist"
assert(dihe_psi), "psi dihedral doesn't exist"


def calc_torsion(dihedral):
    """atom 1 -4 are valid atom selections. torsion in degrees is returned"""
    A = u.select_atoms(dihedral.atom1).positions
    B = u.select_atoms(dihedral.atom2).positions
    C = u.select_atoms(dihedral.atom3).positions
    D = u.select_atoms(dihedral.atom4).positions

    dihe = calc_dihedrals(A, B, C, D)
    return np.rad2deg(dihe)


u = mda.Universe(args.istr, args.itraj,
                 topology_format=args.istrext, format=args.itrajext)

phi_trajdata = np.array(
    [(u.trajectory.frame, calc_torsion(dihe_phi)) for ts in u.trajectory])
psi_trajdata = np.array(
    [(u.trajectory.frame, calc_torsion(dihe_psi)) for ts in u.trajectory])

print(phi_trajdata, psi_trajdata)

phi_frame, phi_series = phi_trajdata.T
psi_frame, psi_series = psi_trajdata.T

phi_series = np.concatenate(phi_series, axis=0)
psi_series = np.concatenate(psi_series, axis=0)

zip(phi_frame, phi_series, psi_series)

with open(args.output, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(phi_frame, phi_series, psi_series))

with sns.axes_style("white"):
    h = sns.jointplot(x=phi_series, y=psi_series, kind="kde", legend=True)
    h.set_axis_labels(r'$\Phi$ (degrees)', r'$\Psi$ (degrees)')
    h.ax_joint.set_xlim(-180, 180)
    h.ax_joint.set_ylim(-180, 180)
    plt.savefig(args.oramachandran_plot, format='png')