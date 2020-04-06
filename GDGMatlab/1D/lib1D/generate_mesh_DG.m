%-------------------------------------------------------------------------
% generate_mesh_DG receives a structure "geo" describing the interval (a,b)
% with mesh size geo.h, and a basis_type variable.
% generate_mesh_DG returns: 
% P (array) containing position of vertices of the mesh,
% T (matrix) where T(:,k) contains position of vertices of element k. 
% P,T do not depend on basis type
% Pb (array) containing position of degrees of freedom,
% Tb (matrix) where Tb(:,k) contains position of degrees of freedom of element k.
% Eb (matrix) where Eb(:,k) contains indexes DOF which are adjacents to
% edge k
% Pb,Eb,Tb do depend on the basis_type.

% author: Tommaso Vanzan
%-------------------------------------------------------------------------

function [P,T,Pb,Tb,Eb]=generate_mesh_DG(geo,basis_type)
NE=(geo.b-geo.a)/geo.h; %Number of elements
%geometrical information
    for k=1:NE+1
        P(k)=geo.a+(k-1)*geo.h; % Matrix with vertices
    end
    for k=1:NE
        T(1,k)=k; % T(:,k) contains vertices of element k.
        T(2,k)=k+1;
    end
%===========
% Linear
%===========
    if  basis_type==101 
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
            
%===========
% Quadratic
%===========        

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