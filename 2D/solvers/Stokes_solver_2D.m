%-------------------------------------------------------------------------
% Stokes_solver_2D receives
% geo: structure containing description of the geometry and
% data/para: structure containing parameters of the problems
% Stokes_solver_2D returns:
% unknowns u,v,p
% Triangle matrix T
% FE matrices Ebve,Pbve,Tbve,Pbpre,Tppre. 

% author: Tommaso Vanzan
%-------------------------------------------------------------------------


function [u,v,p,T,Ebve,Pbve,Tbve,Pbpre,Tbpre]=Stokes_solver_2D(geo,data,para)
%=== Mesh generation
[P,T,E,Pbve,Tbve,Ebve]=generate_mesh_2D(geo,para.basis_type_ve);
Pbpre=P;% DOFs correspond to vertex triangles
Tbpre=T(1:end-1,:);
[Ebpre]=buildmatrixEb_P1_FEM(P,T,E);
%=== Assembly of Stiffness matrix
ndofve=size(Pbve,2);% number of DOF. Supposing conforming FE spaces
ndofpre=size(Pbpre,2);
A1=assemble_matrix_2D(data.mu,P,T,Tbve,Tbve,ndofve,ndofve,para.basis_type_ve,1,0,para.basis_type_ve,1,0,para.order_Gauss);
A2=assemble_matrix_2D(data.mu,P,T,Tbve,Tbve,ndofve,ndofve,para.basis_type_ve,0,1,para.basis_type_ve,0,1,para.order_Gauss);
A3=assemble_matrix_2D(data.mu,P,T,Tbve,Tbve,ndofve,ndofve,para.basis_type_ve,1,0,para.basis_type_ve,0,1,para.order_Gauss);
%A4=para.mu*assemble_matrix_2D(data.mu,P,T,Tbve,Tbve,ndofve,ndofve,para.basis_type_ve,0,1,para.basis_type_ve,1,0,para.order_Gauss);
A5=-assemble_matrix_2D(@(x,y,El) 1,P,T,Tbpre,Tbve,ndofve,ndofpre,para.basis_type_pre,0,0,para.basis_type_ve,1,0,para.order_Gauss);
A6=-assemble_matrix_2D(@(x,y,El) 1,P,T,Tbpre,Tbve,ndofve,ndofpre,para.basis_type_pre,0,0,para.basis_type_ve,0,1,para.order_Gauss);
%A7=-assemble_matrix_2D(@(x,y,El) 1,P,T,Tbve,Tbpre,ndofpre,ndofve,para.basis_type_ve,1,0,para.basis_type_pre,0,0,para.order_Gauss);
%A8=-assemble_matrix_2D(@(x,y,El) 1,P,T,Tbve,Tbpre,ndofpre,ndofve,para.basis_type_ve,0,1,para.basis_type_pre,0,0,para.order_Gauss);
Zero=zeros(ndofpre);
%A=[2*A1+A2,A3,A5;A4,2*A2+A1,A6;A7,A8,Zero];
A=[2*A1+A2,A3,A5;A3',2*A2+A1,A6;A5',A6',Zero];

boundaryu=generate_boundaryedges(Tbve,Ebve,E,data.labelu,para.basis_type_ve);
boundaryv=generate_boundaryedges(Tbve,Ebve,E,data.labelv,para.basis_type_ve);


b=zeros(size(A,1),1);
%=== Assembly RHS
bve1=assemble_rhs_2D(data.f1,P,T,Tbve,ndofve,para.basis_type_ve,0,0); 
bve2=assemble_rhs_2D(data.f2,P,T,Tbve,ndofve,para.basis_type_ve,0,0); 
b(1:2*ndofve)=[bve1;bve2];

%=== Take Care of Boundary

[A,b,flag]=treat_boundary_Stokes(A,b,boundaryu,boundaryv,P,Pbve,Tbve,data,para);

if(flag==0)%if flag==0 I do not have neither Neumann or Robin bc=> need to fix pressure!
    A(2*ndofve+para.index_fix_pressure,:)=0;
    A(2*ndofve+para.index_fix_pressure,2*ndofve+para.index_fix_pressure)=1;
    b(2*ndofve+para.index_fix_pressure)=data.pset(Pbpre(1,para.index_fix_pressure),Pbpre(2,para.index_fix_pressure));
end
 %=== Solution
result=A\b;
u=result(1:ndofve);
v=result(ndofve+1:2*ndofve);
p=result(2*ndofve+1:end);
end