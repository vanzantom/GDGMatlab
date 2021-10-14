function [A,b,w]=treat_boundary_symmetric(A,b,boundary,Pb,c_fun,Dirichlet_fun,Neumann_fun)
% treat_boundary_symmetric takes care of the boundary conditions.
% It eliminates explicitely the Dirichlet nodes, returning a matrix A where
% the row and columns related to a Dirichlet node are eliminated.
% Further it returns an updated b containing the BC, and a vector w which
% contains only the DOF which are still in the reduced A.

%Dirichlet condition: -1
nbn=size(boundary,2);
v=(1:size(A,1));
for k=1:nbn
    if boundary(1,k)==-2 %Neumann case
        i=boundary(2,k);
        b(i)=b(i)+feval(Neumann_fun,Pb(i))*feval(c_fun,Pb(i));
    end
    if boundary(1,k)==-1 %Dirichlet case
        i=boundary(2,k);
        v(i)=-1;
        if(i==1)
            b=b-A*[feval(Dirichlet_fun,Pb(i));zeros(size(A,1)-1,1)];
        else
            b=b-A*[zeros(size(A,1)-1,1);feval(Dirichlet_fun,Pb(i))];
        end
    end
    
end

w=find(v~=-1);
b=b(w);
A=A(w,w);
end