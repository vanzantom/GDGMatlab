function [A,b]=treat_boundary(A,b,boundary,P,Pb,c_fun,Dirichlet_fun,Neumann_fun)
%Dirichlet condition: -1
h=P(2)-P(1);
nbn=size(boundary,2);
for k=1:nbn
    if boundary(1,k)==-2
        i=boundary(2,k);
        b(i)=b(i)+feval(Neumann_fun,Pb(i))*feval(c_fun,Pb(i));
    end
    if boundary(1,k)==-1
        i=boundary(2,k);
        A(i,:)=0;
        A(i,i)=1;
        b(i)=feval(Dirichlet_fun,Pb(i));
    end
    
end
        
    