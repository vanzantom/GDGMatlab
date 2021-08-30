function [ANeumann,A,ARobin,ADir,M_GG,b,boundaryv,e,Pbve,Tbve,Ebve,P,T,E]=generate_matrices_Stokes(geo,data,para)
%generate_matrices_Stokes receives the data structure geo,data and para.
% It returns the matrices
%ANeumann:  stifnesss matrix with zero stress condition everywhere
%A: stiffness matrix with Dirichlet on -1 and Neumann on -2 and -3.
%ADir: stifnesss matrix with Dirichlet BC where imposed on -1 and -3 and
%Neumann on -2 (needed for DtN).
%ARobin: stifness matrix with Dirichlet,Neumann and Robin BC.
%M_GG: local mass matrix on the interior DOF of the interface.
%b: force term which contains Neumann BC and Dirichlet BC.
% e: vector which contains indexes DOfs in interior edge.
% Mass matrix: M.
%Data structures Pbve,Tbve Ebve for unknowns FE and mesh matrix P,T,E.

% author: Tommaso Vanzan
%-------------------------------------------------------------------------
%========================================================================
geo.fun=geo.funStokes;
%=== Mesh generation
[P,T,E,Pbve,Tbve,Ebve]=generate_mesh_2D(geo,para.basis_type_ve);
Pbpre=P;% DOFs correspond to vertex triangles
Tbpre=T(1:end-1,:);
[Ebpre]=buildmatrixEb_P1_FEM(P,T,E);%=== Identification of Boundary
%=== Assembly of Stiffness matrix
ndofve=size(Pbve,2);% number of DOF. Supposing conforming FE spaces
ndofpre=size(Pbpre,2);
A1=assemble_matrix_2D(data.mu,P,T,Tbve,Tbve,ndofve,ndofve,para.basis_type_ve,1,0,para.basis_type_ve,1,0,para.order_Gauss);
A2=assemble_matrix_2D(data.mu,P,T,Tbve,Tbve,ndofve,ndofve,para.basis_type_ve,0,1,para.basis_type_ve,0,1,para.order_Gauss);
A3=assemble_matrix_2D(data.mu,P,T,Tbve,Tbve,ndofve,ndofve,para.basis_type_ve,1,0,para.basis_type_ve,0,1,para.order_Gauss);
%A4=assemble_matrix_2D(data.mu,P,T,Tbve,Tbve,ndofve,ndofve,para.basis_type_ve,0,1,para.basis_type_ve,1,0,para.order_Gauss);
A5=-assemble_matrix_2D(@(x,y,El) 1,P,T,Tbpre,Tbve,ndofve,ndofpre,para.basis_type_pre,0,0,para.basis_type_ve,1,0,para.order_Gauss);
A6=-assemble_matrix_2D(@(x,y,El) 1,P,T,Tbpre,Tbve,ndofve,ndofpre,para.basis_type_pre,0,0,para.basis_type_ve,0,1,para.order_Gauss);
%A7=-assemble_matrix_2D(data.mu,P,T,Tbve,Tbpre,ndofpre,ndofve,para.basis_type_ve,1,0,para.basis_type_pre,0,0,para.order_Gauss);
%A8=-assemble_matrix_2D(data.mu,P,T,Tbve,Tbpre,ndofpre,ndofve,para.basis_type_ve,0,1,para.basis_type_pre,0,0,para.order_Gauss);
Zero=zeros(ndofpre);
%ANeumann=[2*A1+A2,A3,A5;A4,2*A2+A1,A6;A7,A8,Zero];
ANeumann=[2*A1+A2,A3,A5;A3',2*A2+A1,A6;A5',A6',Zero];

boundaryu=generate_boundaryedges(Tbve,Ebve,E,data.labelu,para.basis_type_ve);
boundaryv=generate_boundaryedges(Tbve,Ebve,E,data.labelv,para.basis_type_ve);
b=zeros(size(ANeumann,1),1);
%=== Assembly RHS
bve1=assemble_rhs_2D(data.f1,P,T,Tbve,ndofve,para.basis_type_ve,0,0); 
bve2=assemble_rhs_2D(data.f2,P,T,Tbve,ndofve,para.basis_type_ve,0,0); 
b(1:2*ndofve)=[bve1;bve2];
%=== Take Care of Boundary
[A,ARobin,ADir,b,e,M_GG]=impose_boundary_Stokes(ANeumann,b,boundaryu,boundaryv,P,Pbve,Tbve,Tbve,data,para);
if(para.index_fix_pressure>0)% if para.index_fix_pressure>0, then I have Dirichlet everywhere and I need to fix pressure
    ADir(2*ndofve+para.index_fix_pressure,:)=0;
    ADir(2*ndofve+para.index_fix_pressure,2*ndofve+para.index_fix_pressure)=1;
end
end