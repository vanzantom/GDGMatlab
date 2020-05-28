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



% author: Tommaso Vanzan
%-------------------------------------------------------------------------

function [A,b,flag]=treat_boundary_Stokes(A,b,boundaryu,boundaryv,P,Pb,Tb_test,data,para)

%========================================================================
nbn=size(boundaryu,2);
ndofve=size(Pb,2);
Tb_trial=Tb_test;
number_of_local_basis_test=size(Tb_test,1); %number of local basis function on element El
number_of_local_basis_trial=size(Tb_trial,1); %number of local basis function on element El
flag=0;
for k=1:nbn
    if boundaryv(1,k)==-2
            flag=1;
            El=boundaryu(2,k);% Element which as edge boundary(2,k).
            number_of_local_basis_test=size(Tb_test,1); %number of local basis function on element El
            Vertex_El=P(:,Tb_test(1:3,El));%vertices of the element
            vertices=P(:,[boundaryv(3,k);boundaryv(4,k)]);%vertices of the edge
            edge_length=norm(vertices(:,1)-vertices(:,2),2);
            normal=find_normal(boundaryv(3,k),boundaryv(4,k),Tb_test(1:3,El),P);
            [~,i]=max(normal);
            sign=(-1)^i;
            [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(vertices(:,1),vertices(:,2),para.order_Gauss);%already multiplied by length edge
                for alpha=1:number_of_local_basis_test
                   b(Tb_test(alpha,El)+ndofve)=b(Tb_test(alpha,El)+ndofve)+ sign*Gauss_quadrature_1D_test(data.Normal_stress,Gauss_nodes',Gauss_weights,Vertex_El,para.basis_type_ve,alpha,0,0);
                end
    elseif boundaryv(1,k)==-3 %Robin Boundary condition!
            flag=1;
            El=boundaryv(2,k);% Element which as edge boundary(2,k).
            number_of_local_basis_test=size(Tb_test,1); %number of local basis function on element El
            Vertex_El=P(:,Tb_test(1:3,El));%vertices of the element
            vertices=P(:,[boundaryv(3,k);boundaryv(4,k)]);%vertices of the edge
            edge_length=norm(vertices(:,1)-vertices(:,2),2);
            normal=find_normal(boundaryv(3,k),boundaryv(4,k),Tb_test(1:3,El),P);
            [~,i]=max(normal);
            sign=(-1)^i;
            [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(vertices(:,1),vertices(:,2),para.order_Gauss);%already multiplied by length edge
                for alpha=1:number_of_local_basis_trial %loop over trial functions on the element for trial FE space.
                    for beta=1:number_of_local_basis_test %loop over trial functions on the element for test FE space.
                        int_value=Gauss_quadrature_1D_trial_test(@(x,y) para.Robin,Gauss_nodes,Gauss_weights,Vertex_El,Vertex_El,para.basis_type_ve,alpha,0,0,para.basis_type_ve,beta,0,0);
                        A(Tb_test(beta,El)+ndofve,Tb_trial(alpha,El)+ndofve)=A(Tb_test(beta,El)+ndofve,Tb_trial(alpha,El)+ndofve)+para.Robin*int_value;
                    end
                    b(Tb_test(alpha,El)+ndofve)=b(Tb_test(alpha,El)+ndofve)+ sign*Gauss_quadrature_1D_test(data.Robin_stress,Gauss_nodes',Gauss_weights,Vertex_El,para.basis_type_ve,alpha,0,0);
                end
    end
    
end
for k=1:nbn
    if boundaryu(1,k)==-1
       index=boundaryu(3:end,k);% indexes of nodes which are on the edge boundary
        for j=1:length(index)
            A(index(j),:)=0;
            A(index(j),index(j))=1;
            b(index(j))=feval(data.Dirichlet_u,Pb(1,index(j)),Pb(2,index(j)));
        end
    end
     if boundaryv(1,k)==-1
       index=boundaryv(3:end,k);% indexes of nodes which are on the edge boundary
       for j=1:length(index)
            A(index(j)+ndofve,:)=0;
            A(index(j)+ndofve,index(j)+ndofve)=1;
            b(index(j)+ndofve)=feval(data.Dirichlet_v,Pb(1,index(j)),Pb(2,index(j)));
       end
     end
end
end
    