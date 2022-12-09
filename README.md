# Computational Molecular Dynamics

A simple computational code to compute the basic behaviour of a system of particles and determine its consistency under the Lennard Jones Potential Curve. 

## Description
One of the most widely used intermolecular potentials in classical many-body simulations, is the so-called Lennard-Jones 12-6 potential,
> $\displaystyle V(r) = 4 \epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]$

Also expressed as,
> $\displaystyle V(r) = \left(\frac{A}{r^{12}}\right) - \left(\frac{B}{r^{6}}\right)$\

where,

1. $V$ is the intermolecular potential between the two atoms or molecules.
2. $\epsilon$ is the well depth and a measure of how strongly the two particles attract each other.
3. $\sigma$ is the distance at which the intermolecular potential between the two particles is zero (Figure 1). $\sigma$ gives a measurement of how close two non-bonding particles can get and is thus referred to as the van der Waals radius. It is equal to one-half of the internuclear distance between nonbonding particles.
4. $r$ is the distance of separation between both particles (measured from the center of one particle to the center of the other particle).
5. $A = 4 \epsilon \sigma^{12}$, $B = 4 \epsilon \sigma^{6}$ 
6. Minimum value of $\Phi(r)$ at $r=r_{min}$.

Lennard-Jones-type $r{−n}-r^{−m}$ pair potentials were proposed in 1925 by Jones [[1]](<https://royalsocietypublishing.org/doi/10.1098/rspa.1924.0081> "J. Jones, Proc R Soc London A, 1924, 106, 441–462") (later Lennard-Jones) to describe the cohesive energy of crystals of noble gases, such as Argon. The now conventional LJ 12-6 form was proposed by Lennard-Jones in 1931 [[2]](<https://iopscience.iop.org/article/10.1088/0959-5309/43/5/301> "J. E. Lennard-Jones, Proc Phys Soc, 1931, 43, 461–482") after London had derived that the dispersion interaction between atoms decays as $r^{−6}$ (at least, in the non-retarded regime). [[3]](<https://arxiv.org/pdf/1910.05746.pdf> "Wang, Xipeng, et al. "The Lennard-Jones potential: when (not) to use it." Physical Chemistry Chemical Physics 22.19 (2020): 10624-10633")

## Getting Started

### Dependencies

The code exists in two variants, C++ and Python and do not require additional installation other than the basic language compilers.

### Executing program

* C++
```
cd /Home/User/File_Location
g++ Molecular_Dynamics.cpp
```

* Python
```
cd /Home/User/File_Location
python3 Translated_MD.cpp
```
