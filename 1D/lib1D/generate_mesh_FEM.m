%-------------------------------------------------------------------------
% generate_mesh_FEM receives a structure "geo" describing the interval (a,b)
% with mesh size geo.h, and a basis_type variable.
% generate_mesh_FEM returns: 
% P (array) containing position of vertices of the mesh,
% T (matrix) where T(:,k) contains position of vertices of element k.
% Pb (array) containing position of degrees of Freedom elements,
% Tb (matrix) where Tb(:,k) contains position of degrees of freedom of element k.

% author: Tommaso Vanzan
%-------------------------------------------------------------------------

function [P,T,Pb,Tb]=generate_mesh_FEM(geo,basis_type)
NE=(geo.b-geo.a)/geo.h; %Number of elements
%Mesh information
    for k=1:NE+1
        P(k)=geo.a+(k-1)*geo.h; % Matrix with vertices
    end
    for k=1:NE
        T(1,k)=k; % T(:,k) contains vertices of the k elements
        T(2,k)=k+1;
    end
% Create DOF matrices according to basis_type
    if basis_type==101 % Linear FE
        Pb=P;
        Tb=T;
    elseif basis_type==102 %Quadratic FE
        for k=1:2*NE+1
            Pb(k)=geo.a+(k-1)*geo.h/2; % Matrix with vertices
        end
        j=1;
        for k=1:NE
            Tb(1,k)=j; % T(:,k) contains vertices of the k elements
            Tb(2,k)=j+1;
            Tb(3,k)=j+2;
            j=j+2;
        end
    end
end