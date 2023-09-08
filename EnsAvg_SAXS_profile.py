# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 11:35:05 2022

@author: popara
"""

import sys
sys.path.insert(1,'C:\\Program Files\\IMP-2.17.0\\python\\') 
import IMP
import IMP.atom
import IMP.core
import IMP.saxs


import numpy as np
import pathlib
import os
import tqdm

import matplotlib as mpl
import matplotlib.pyplot as plt




def get_model_profile(
    pdb_fn: pathlib.Path, 
    model_delta_q = 0.5 / 500,
    model_min_q = 0.0,
    model_max_q = 0.5
) -> IMP.saxs.Profile:
    m = IMP.Model()
    pdb_fn = str(pdb_fn)
    saxs_fn = pdb_fn + '.saxs.dat'
    # calculate SAXS profile
    model_profile = IMP.saxs.Profile(model_min_q, model_max_q, model_delta_q)
    if os.path.exists(saxs_fn):
        model_profile.read_SAXS_file(saxs_fn)
    else:
        mp = IMP.atom.read_pdb(pdb_fn, m, IMP.atom.NonWaterNonHydrogenPDBSelector(), True, True)
        # select particles from the model
        particles = IMP.atom.get_by_type(mp, IMP.atom.ATOM_TYPE)
        # add radius for water layer computation
        ft = IMP.saxs.get_default_form_factor_table()
        for i in range(0, len(particles)):
            radius = ft.get_radius(particles[i])
            IMP.core.XYZR.setup_particle(particles[i], radius)
        # compute surface accessibility
        s = IMP.saxs.SolventAccessibleSurface()
        surface_area = s.get_solvent_accessibility(IMP.core.XYZRs(particles))
        model_profile.calculate_profile_partial(particles, surface_area)
        # Write SAXS curve
        model_profile.write_SAXS_file(saxs_fn)
    return model_profile




######## input parameters and data paths ##########################



pdb_path = pathlib.Path('C:/user/SAXS_profile/example_data/PDBs/')

weights_path = pathlib.Path('C:/user/SAXS_profile/example_data/')

outfile_path = pathlib.Path('C:/user/SAXS_profile/')


# q vector range for calculation of theoretical scattering profiles
model_delta_q = 0.6 / 500 
model_min_q = 0.0
model_max_q = 0.6  



################### Compute SAXS profiles for PDBs and ensemble average ##########



weights = np.loadtxt(weights_path/'conformer_weights.dat', delimiter=' ', usecols=[1])

pdb_fns = sorted(list(pdb_path.glob('*.pdb')))

assert len(weights)==len(pdb_fns), "Length of weights file does not match number of conformers"



model_profiles = list()
for fn in tqdm.tqdm(pdb_fns):
    model_profiles.append(
        get_model_profile(fn, model_delta_q, model_min_q, model_max_q)
    )

   

avg_saxs_fn = str(outfile_path / 'ensemble_averaged_saxs_profile.dat')

EnsAvg_profile = IMP.saxs.Profile(model_min_q, model_max_q, model_delta_q)
for p, w in zip(model_profiles, weights):
    EnsAvg_profile.add(p, w)


EnsAvg_profile.write_SAXS_file(avg_saxs_fn)



################ Visualization #################################



mpl.rcParams['font.sans-serif']="Arial"
plt.rcParams['figure.figsize'] = 3.25, 4.5
params = {'legend.fontsize': 6,
          'legend.handlelength': 2}

plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fancybox'] = True
mpl.rcParams['axes.linewidth'] = 0.3

mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['xtick.labelsize'] = 6
mpl.rcParams['ytick.labelsize'] = 6
mpl.rcParams['xtick.major.width'] = 0.3
mpl.rcParams['ytick.major.width'] = 0.3
mpl.rcParams['xtick.minor.width'] = 0.3
mpl.rcParams['ytick.minor.width'] = 0.3
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['xtick.minor.size'] = 1 # half of the major ticks length
mpl.rcParams['ytick.minor.size'] = 1

plt.rcParams.update(params)




fig = plt.figure() 
gs = fig.add_gridspec(2, 1, hspace=0, wspace=0) 
ax_t = fig.add_subplot(gs[0,0]) # top panel I versus q
ax_b = fig.add_subplot(gs[1,0]) # bottom; Kratky plot I*q^2 versus q 

    

saxs_data = np.loadtxt(avg_saxs_fn)
    
ax_t.plot(saxs_data[:,0],saxs_data[:,1],'-', lw=0.8)
ax_b.plot(saxs_data[:,0],saxs_data[:,1]*saxs_data[:,0]**2,'-', lw=0.8)
    
        
ax_t.set_yscale('log')
ax_t.set_ylabel(r'$I(q)$')
ax_t.set_xlim(model_min_q, None)
ax_b.set_xlim(model_min_q, None)
ax_t.axes.xaxis.set_ticklabels([])
ax_b.set_xlabel(r'$q$ $[Å^{-1}$]')
ax_b.set_ylabel(r'$I(q) \cdot q^{2}$')



fig.tight_layout()
fig.savefig(outfile_path/'EnsAvg_saxs_profile.png',dpi=300, transparent=True)
fig.savefig(outfile_path/'EnsAvg_saxs_profile.svg', transparent=True)
    
plt.close(fig)

    