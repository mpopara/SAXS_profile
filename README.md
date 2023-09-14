# SAXS_profile
Python script to compute theoretical ensemble-averaged SAXS profile for a biomolecule


## General description

This script computes theoretical Small Angle X-ray Scattering (SAXS) profiles for the ensemble of conformational models, as well as the ensemble-averaged SAXS profile.<sup>1</sup> 
SAXS profiles are computed using _saxs_ module<sup>2</sup> of the Integrative Modelling Platform<sup>3</sup> (IMP), for a user-defined range of momentum transfer vector _q_ values, from _q_<sub>min</sub> to _q_<sub>max</sub>, with a step of _dq_.
Ensemble-averaged SAXS profile is computed as weighted linear combination of SAXS profiles of individual ensemble members. As weights used are population fractions of conformational models in the ensemble.

Individual and ensemble-averaged SAXS profiles are saved as .dat file, and furthermore, results are visualized and saved as png and svg file.
Figure comprises of two subplots (see below). Top subplot contains the scattering profile, i.e. scattering intensity _I_(_q_) verus the momentum transfer vector _q_. 
Shape of the scattering profile informs about the flexibility of the biomolecule<sup>4</sup>, i.e. folded proteins have scattering profile with rises and dips, while unfolded proteins have featureless scattering profile.
In between those two scenarios are multi-domain proteins, which have profile with smoothed features. 
Lover panel contains I(_q_)*_q_<sup>2</sup> versus _q_, known as Kratky plot. Its shape informs on the compactness of the biomolecule, and goes from parabolic shape for folded proteins, to hyperbolic for unfolded proteins.<sup>4</sup>


![EnsAvg_saxs_profile](https://github.com/mpopara/SAXS_profile/assets/40856779/a05738f4-d4a9-4b44-8b73-4c522a63018b)


## Example data and input file requirements

In the folder example_data, provided are exemplary input files:

* /PDBs/*.pdb- ensemble of structural models, provided as individual .pdb files.
* .dat file containing weights (population fractions) of ensemble members. This space-delimited file is of a size N<sub>conformers</sub> x 2, where the first column contains indices of the ensemble members,
 and the second column contains their corresponding weights. This script assumes that the numbering of esemble members in their file name follows the same order as in the weights file. 

## Dependencies

EnsAvg_SAXS_profile.py is a python script bult on Python 3.8.8. Script was tested with provided examplary input files under the following configuration:

* Windows 10
* Python 3.8.8
* IMP 2.17.0
* numpy 1.23.0
* matplotlib 3.7.1

## References
1. Dittrich, J.; Popara, M.; Kubiak, J.; Dimura, M.; Schepers, B.; Verma, N.; Schmitz,
B.; Dollinger, P.; Kovacic, F.; Jaeger, K. E.; Seidel, C. A. M.; Peulen, T. O.; Gohlke, H.,
Resolution of Maximum Entropy Method-Derived Posterior Conformational Ensembles of a
Flexible System Probed by FRET and Molecular Dynamics Simulations. J Chem Theory
Comput 2023, 19 (8), 2389-2409.

2. Schneidman-Duhovny, D.; Hammel, M.; Tainer, John A.; Sali, A., Accurate SAXS Profile Computation and its Assessment by Contrast Variation Experiments. Biophys. J. 2013, 105 (4), 962-974.

3. Russel, D.; Lasker, K.; Webb, B.; Velázquez-Muriel, J.; Tjioe, E.; Schneidman-Duhovny, D.; Peterson, B.; Sali, A., Putting the Pieces Together: Integrative Modeling Platform Software for Structure Determination of Macromolecular Assemblies. PLoS Biol. 2012, 10 (1), e1001244.

4. Brosey, C. A.; Tainer, J. A., Evolving SAXS versatility: solution X-ray scattering for macromolecular architecture, functional landscapes, and integrative structural biology. Curr. Opin. Struct. Biol. 2019, 58, 197-213.


## Authors

* Thomas-Otavio Peulen
* Milana Popara
