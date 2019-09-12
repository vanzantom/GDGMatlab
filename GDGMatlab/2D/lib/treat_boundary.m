function [A,b]=treat_boundary(A,b,boundary,Pb,Dirichlet_fun)
%treat_boundary is a function which receives the Matrix A, the right hand
%side b, the matrix of the DOF Pb, and the Dirichlet data Dirichlet_fun.
%According to the label associated to the boundary nodes, it implements
%Dirichlet BC(-1) or Neumann(-2)

% NB: Neumann BC are still to implement!
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
        
    