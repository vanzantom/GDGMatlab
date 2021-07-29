%-------------------------------------------------------------------------
% Stokes_Darcy_solver_2D receives
% geo: structure containing description of the geometry and
% data/para: structure containing parameters of the problems both for
% Stokes and Darcy

% Stokes_solver_2D returns:
% results containing solution
% Atot: global matrix
% C: coupling matrix between Stokes and Darcy
% matrixS and matrixD: structure containing matrices Stokes and Darcy
% corresponding to different BC.
% FES/FED structures referring FE info for both Stokes and Darcy


% author: Tommaso Vanzan
%-------------------------------------------------------------------------


function [result,Atot,C,matrixS,FES,matrixD,FED]=Stokes_Darcy_solver_2D(geo,dataStokes,paraStokes,dataDarcy,paraDarcy)
%Build Global matrix
geo.fun=geo.funDarcy;
[matrixD.ANeumann,matrixD.A,matrixD.ARobin,matrixD.ADir,matrixD.M_GG,bD,FED.eD,FED.Pb,FED.Tb,FED.P,FED.T,FED.E,FED.boundary]=generate_matrices_Darcy(geo,dataDarcy,paraDarcy);
geo.fun=geo.funStokes;
[matrixS.ANeumann,matrixS.A,matrixS.ARobin,matrixS.ADir,matrixS.M_GG,bS,FES.boundaryv,FES.eS,FES.Pb,FES.Tb,FES.Eb,FES.P,FES.T,FES.E]=generate_matrices_Stokes(geo,dataStokes,paraStokes);
ndofve=size(FES.Pb,2);%# dof velocities
ndofpre=size(FES.P,2);%# dof Stokes pressure
ndofq=size(FED.Pb,2);%# dof Darcy pressyre
[C,C1,C2]=assembly_coupling_matrix(FES.Pb,FES.Tb,FED.Pb,FED.Tb,FES.boundaryv,FED.boundary,geo.label,paraStokes,paraDarcy,0,1,0);
Atot=blkdiag(matrixS.A,matrixD.A); %assemble coupled system
Atot(ndofve+1:2*ndofve,2*ndofve+ndofpre+1:end)=C1;%Insert coupling terms
Atot(2*ndofve+ndofpre+1:end,ndofve+1:2*ndofve)=C2;


%=== Solve
result=Atot\[bS;bD];


end
