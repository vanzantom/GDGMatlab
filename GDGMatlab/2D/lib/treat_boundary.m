%-------------------------------------------------------------------------
% treat_boundary receives 
% A: stiffness matrix
% b: right hand side vector
% boundary: matrix describing the boundary conditions on the boundary
% Pb: matrix with position degrees of freedom.
% Dirichlet_fun: Dirichlet boundary data
%treat_boundary returns
% A: modified stiffness matrix where the Dirichlet BC are strongly imposed
% b: modified right hand side vector to include boundary conditions

% Attention: Neumann BC are still to implement!

% author: Tommaso Vanzan
%-------------------------------------------------------------------------

function [A,b]=treat_boundary(A,b,boundary,Pb,Dirichlet_fun)

%========================================================================
nbn=size(boundary,2);
for k=1:nbn
    if boundary(1,k)==-1
        i=boundary(2,k);
        A(i,:)=0;
        A(i,i)=1;
        b(i)=feval(Dirichlet_fun,Pb(1,i),Pb(2,i));
    end
    
end
        
    