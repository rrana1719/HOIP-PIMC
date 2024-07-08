# HOIP-PIMC

This repository has source code and figure code for "On the interplay of electronic and lattice screening on exciton binding in two-dimensional lead halide perovskites". We aim to understand exciton-phonon interactions in 2D hybrid organic-inorganic perovskites (HOIPs).

We looked into two geometries.

1) The folder titled "Crystal Geometry" has the electron-hole potential (Hansen et al.) for an exciton in an infinitely repeating lattice of organic and inorganic layers (titled "n1.txt"). The python script in the folder uses path integral Monte Carlo (PIMC) with a Fröhlich-like model Hamiltonian to get the renormalized exciton binding energy and corresponding P(r) at a certain temperature. We note that the code can get the exciton binding energy in the absence of phonon interactions if needed.
   
2) The folder titled "Thin Film Geometry" has everything incorporated into the Python script. The Hankel Transform library is used to get the electron-hole potential for an exciton in a thin inorganic film surrounded by semi-infinite organic layers. Like the crystal geometry, the python script uses path integral Monte Carlo (PIMC) with a Fröhlich-like model Hamiltonian to get the renormalized exciton binding energy and corresponding P(r) at a certain temperature. We note that this code can also get the exciton binding energy in the absence of phonon interactions if needed.
   


Credit to:
K. R. Hansen, C. Y. Wong, C. E. McClure, B. Romrell, L. Flannery, D. Powell, K. Garden, A. Berzansky, M. Eggleston, D. J. King,
et al., Mechanistic origins of excitonic properties in 2d perovskites: Implications for exciton engineering, Matter 6, 3463 (2023).
