#!/usr/bin/env python

import argparse
import csv
import sys

import itertools
import MDAnalysis as mda

import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt

import numpy as np
import numpy.linalg

def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--itraj', help='input traj')
    parser.add_argument('--istr', help='input str')
    parser.add_argument('--itrajext', help='input traj ext')
    parser.add_argument('--istrext', help='input str ext')
    parser.add_argument('--isegid1', help='segid 1')
    parser.add_argument('--ilabel', help='plot label')
    parser.add_argument('--ititle1', help='plot title')
    parser.add_argument('--output1', help='output1 - timeseries')
    parser.add_argument('--o_plot', help='End to End plot')
    return parser.parse_args()


args = parse_command_line(sys.argv)



u = mda.Universe(args.istr, args.itraj,
                 topology_format=args.istrext, format=args.itrajext)

bbatoms = "(segid %s and backbone)" % \
    (args.isegid1)

bb = u.select_atoms(bbatoms)  # a selection (a AtomGroup)
bbfirst = u.select_atoms(bbatoms)[0]
bblast = u.select_atoms(bbatoms)[-1]



rgyr = []

for ts in u.trajectory:  # iterate through all frames
    rg = bb.radius_of_gyration()  # method of a AtomGroup; updates with each frame
    rgyr.append((ts.frame, rg))

rgyr = np.array(rgyr)


color = itertools.cycle(['r', 'b', 'gold'])

fig, axs = plt.subplots(1, 2, sharex=False, sharey=False, tight_layout=True)

params = {
   'axes.labelsize': 8,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'figure.figsize': [4.5, 4.5],
   'figure.dpi':300
   }
plt.rcParams.update(params)

axs[0].plot(rgyr[:,0], rgyr[:,1], 'r-', lw=2, label=args.ilabel)
axs[0].set_xlabel("number of frames")
axs[0].set_ylabel(r"radius of gyration  ($\AA$)")
axs[0].legend()

n, bins, patches = axs[1].hist(rgyr[:,1], color=next(color), label=args.ilabel, alpha=0.5, density=True, stacked=True) #, bins=20)

axs[1].legend()
axs[1].set_ylabel('Density Normalised Frequency');
axs[1].set_xlabel(r'radius of gyration ($\AA$)')
fig.suptitle(args.ititle1, fontsize=12, fontweight='bold')
fig.subplots_adjust(top=0.45)

print(" \n".join(['The Radius of gyration for the backbone atoms of the segid %s. The first and last atoms of the selection are %s %s:',str(args.isegid1),str(bbfirst), str(bblast)]))

plt.savefig(args.o_plot, format='png') # svg is better but sticking with png for now


np.savetxt(args.output1, rgyr, delimiter='\t')
