# GDGMatlab
GDGMatlab contains libraries to implement Finite element methods (FEM) and Discontinous Galerking finite element methods (DG) for the solution of partial differential equations.

PREMISE: The library is not oriented to high performance computing. Througout the codes, a modular, easy to undestand and to generalize approach is preferred over the performance. Thus, GDGMatlab is oriented to educational and research purposes only. 


GDGMatlab is inspired by several works:

- GDG90: A Fortran library, developed by Soheil Hajian, which contains some of the features of GDGMatlab and further ones.
- Lectures notes on the implementation of Finite element methods by Prof. Xiaoming He, Missouri University of Science & Technology.
- Lectures notes by Prof. Antonietti, Politecnico di Milano, on the implementation of DG methods.

 
All possible bugs/errors are my full responsibiliy. 
If you find one, please send an email to tommaso.vanzan@unige.ch.

Currently implemented are:
1D  
-FEM (P1,P2) (Dirichlet and Neumann Boundary conditions (BC))  
-DG-Symmetric Interior Penalty method (P1,P2) (IP)  (homogeneous Dirichlet BC)  
-DG-Hybridizable Symmetric Interior Penalty method (IPH) (P1,P2)  (homogeneous Dirichlet BC)  

2D  
-FEM (P1,P1b,P2,P2b) (Dirichlet/Neumann/Robin BC).  
-P1-DG Symmetric Interior Penalty method (IP). (Homogeneous Dirichlet BC)  
-P1-DG Hybridizable Symmetric Interior Penaly method(IPH). (Homogeneous Dirichlet BC)  
-P1-HDG Hybrid Discontinous Galerking method (HDG). (Homogeneous Dirichlet BC)  

Further it includes example files to solve Poisson,Stokes and Stokes-Darcy problems.


Future developments:

-Implementation of general Dirichlet and Neumann BC in DG 1D.  
-Implementation of general Dirichlet BC in 2D for DG.  
