# GDGMatlab
GDGMatlab contains libraries to implement Finite element methods (FEM) and Discontinous Galerking finite element methods (DG) for the solution of a Poisson Problem.

PREMISE: The library is not oriented to high performance computing. Througout the codes, a modular,easy to undestand and to generalize approach is preferred over the performance. Thus, GDGMatlab is oriented to educational and research purposes only. 


GDGMatlab is inspired by several works. I would like to thank:

-GDG90: A Fortran library, developed by Soheil Hajian, which contains some of the features of GDGMatlab and further ones.
-Lectures notes on the implementation of Finite element methods by Prof. Xiaoming He, Missouri University of Science & Technology.
-Lectures notes by Prof. Antonietti, Politecnico di Milano, on the implementation of DG methods.

 
All possible bugs/errors are my full responsibiliy. 
If you find one, please send an email to tommaso.vanzan@unige.ch.

Currently implemented are:

1D

-P1-FEM (With Dirichlet and Neumann Boundary conditions (BC))

-P2-FEM (With Dirichlet and Neumann BC)

-P1-DG  (With homogeneous Dirichlet BC)

-P2-DG  (With homogeneous Dirichlet BC)

2D

-P1-FEM (General Dirichlet BC).

-P1DG Symmetric Interior Penalty method (IP). (Homogeneous Dirichlet BC)

-P1DG Hybridizable Symmetric Interior Penaly method(IPH). (Homogeneous Dirichlet BC)

-P1HDG Hybrid Discontinous Galerking method (HDG). (Homogeneous Dirichlet BC)

-A further code to implement Optimized Schwarz methods for Hybrid discontinous Galerkin discretization.

Martin J. Gander and Soheil Hajian, Analysis of Schwarz methods for a hybridizable discontinuous Galerkin discretization, SIAM Journal on Numerical Analysis, 53 (2015), pp. 573-597

Future developments:

-Implementation of general Dirichlet and Neumann BC in DG 1D.

-Implementation of Neumann BC in 2D for FEM.

-Implementation of general Dirichlet BC in 2D for DG.
