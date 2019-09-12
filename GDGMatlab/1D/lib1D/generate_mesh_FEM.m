function [P,T,Pb,Tb]=generate_mesh_FEM(a,b,h,basis_type)
% generate_mesh_FEM creates the matrices P,T with mesh infos and matrices Pb,Tb with
% the information of the degree of freedoms.
NE=(b-a)/h; %Number of elements
%Mesh information
    for k=1:NE+1
        P(k)=a+(k-1)*h; % Matrix with vertices
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
            Pb(k)=a+(k-1)*h/2; % Matrix with vertices
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