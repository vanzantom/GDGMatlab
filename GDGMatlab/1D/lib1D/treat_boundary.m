%-------------------------------------------------------------------------
% treat_boundary receives 
% A: stiffness matrix
% b: right hand side vector
% boundary: matrix describing the boundary conditions on the boundary
% vertices
% Pb: matrix with position degrees of freedom.
% c_fun: diffusion coefficient
% Dirichlet_fun,Neumann_fun: Dirichlet/Neumann boundary data
%treat_boundary returns
% A: modified stifness matrix where the Dirichlet BC are strongly imposed
% b: modified right hand side vector to include boundary conditions

% author: Tommaso Vanzan
%-------------------------------------------------------------------------
function [A,b]=treat_boundary(A,b,boundary,P,Pb,c_fun,Dirichlet_fun,Neumann_fun)
%Dirichlet condition: -1
h=P(2)-P(1);%mesh size
nbn=size(boundary,2);
for k=1:nbn
    if boundary(1,k)==-2 %Neumann case
        i=boundary(2,k);
        b(i)=b(i)+feval(Neumann_fun,Pb(i))*feval(c_fun,Pb(i));
    end
    if boundary(1,k)==-1 %Dirichlet case
        i=boundary(2,k);
        A(i,:)=0;
        A(i,i)=1;
        b(i)=feval(Dirichlet_fun,Pb(i));
    end
    
end
        
    