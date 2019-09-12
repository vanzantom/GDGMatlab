function [P,T,Pb,Tb,Eb]=generate_mesh_DG(a,b,h,basis_type)
% generate_mesh_DG creates the matrices P,T with mesh infos and matrices Pb,Tb with
% the information of the degree of freedoms for the DG discretization.
NE=(b-a)/h; %Number of elements
%geometrical information
    for k=1:NE+1
        P(k)=a+(k-1)*h; % Matrix with vertices
    end
    for k=1:NE
        T(1,k)=k; % T(:,k) contains vertices of the k elements
        T(2,k)=k+1;
    end
% Create DOF matrices according to basis_type
    if  basis_type==101 % Linear FE
        n=size(P,2);
        Pb(1)=P(1);
        j=2;
        for k=2:n-1   %two DOFs for each interior node.
            Pb(j)=P(k);
            Pb(j+1)=P(k);
            j=j+2;
        end
        Pb(n*2-2)=P(n);
        j=1;
        for k=1:NE    
            Tb(1,k)=j; % T(:,k) contains vertices of the k elements
            Tb(2,k)=j+1;
            j=j+2;
        end 
        j=1;
        for k=1:NE-1 % I create a matrix of interfaces 'edges'. Contains indices T that share edges
            Eb(1,k)=j;
            Eb(2,k)=j+1;
            j=j+1;
        end
            
        
        %To fix the quadratic DG
    elseif basis_type==102 %Quadratic FE
      n=size(P,2);
        
        j=1;
        for k=1:n-1 %two DOFs for each boundary node of a element plus the middle point.
            Pb(j)=P(k);
            Pb(j+1)=(P(k)+P(k+1))/2;
            Pb(j+2)=P(k+1);
            j=j+3;
        end
        Pb(n+n-2+n-1)=P(n);
        j=1;
        for k=1:NE
            Tb(1,k)=j; % T(:,k) contains vertices of the k elements
            Tb(2,k)=j+1;
            Tb(3,k)=j+2;
            j=j+3;
        end 
        j=1;
        for k=1:NE-1 % I create a matrix of interfaces 'edges'. Contains T that share edges
            Eb(1,k)=j;
            Eb(2,k)=j+1;
            j=j+1;
        end
    end
end