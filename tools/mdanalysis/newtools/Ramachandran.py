import MDAnalysis 
from MDAnalysis.tests.datafiles import PSF,DCD,GRO,XTC  # test trajectory
import numpy.linalg
import numpy as np

from MDAnalysis.analysis.dihedrals import Dihedral
from MDAnalysis.analysis.dihedrals import Ramachandran

#u = MDAnalysis.Universe("step3_charmm2omm.psf", "step5_1.dcd")  
u = MDAnalysis.Universe(PSF,DCD)
r = u.select_atoms("resid 5-10")
#r = u.select_atoms("segid PROA")

R = Ramachandran(r).run()

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')  # noqa
#fig, ax = plt.subplots(figsize=plt.figaspect(1))
#R.plot(ax=ax, color='k', marker='s')
#plt.savefig("rtest.png", format='png') # svg is better but sticking with png for now
#plt.show()


u = MDAnalysis.Universe(GRO, XTC)
r = u.select_atoms("resid 5-10")

R = Ramachandran(r).run()

fig, ax = plt.subplots(figsize=plt.figaspect(1))
R.plot(ax=ax, color='k', marker='s')
plt.savefig("rtest.png", format='png') # svg is better but sticking with png for now

#with sns.axes_style("white"):
#    h = sns.jointplot(x=phi_series, y=psi_series, kind="kde", legend=True)
#    h.set_axis_labels(r'$\Phi$ (degrees)', r'$\Psi$ (degrees)')
#    h.ax_joint.set_xlim(-180, 180)
#    h.ax_joint.set_ylim(-180, 180)
#    plt.savefig(args.oramachandran_plot, format='png')

