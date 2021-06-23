function [A,b,w]=treat_boundary_symmetric(A,b,boundary,Pb,Tb_trial,Tb_test,para,Dirichlet_fun,Neumann_fun,Robin_fun)
%-------------------------------------------------------------------------
% treat_boundary_symmetric receives 
% A: stiffness matrix
% b: right hand side vector
% boundary: boundaryedge matrix
% Pb: matrix with position degrees of freedom.
% Tb_trial: matrix with elements information for trial FE space
% Tb_test: matrix with elements information for test FE space
% Dirichlet_fun: Dirichlet boundary data
% Neumann_fun: Neumann boundary data
% Robin_fun: Neumann boundary data
%treat_boundary returns
% A: modified stiffness matrix where the Dirichlet nodes are eliminated.
% b: modified right hand side vector to include boundary conditions

% author: Tommaso Vanzan
%-------------------------------------------------------------------------
%========================================================================
nbn=size(boundary,2);
v=(1:size(A,1))';
for k=1:nbn
    if boundary(1,k)==-2 %Neumann data
        El=boundary(2,k);% Element which as edge boundary(2,k).
        number_of_local_basis_test=size(Tb_test,1); %number of local basis function on element El
        Vertex_El=Pb(:,Tb_test(:,El));%vertices of the element
        vertices=Pb(:,[boundary(3,k);boundary(4,k)]);%vertices of the edge
        edge_length=norm(vertices(:,1)-vertices(:,2),2);
        [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(vertices(:,1),vertices(:,2),para.order_Gauss);%already multiplied by length edge
        for alpha=1:number_of_local_basis_test
           b(Tb_test(alpha,El))=b(Tb_test(alpha,El))+ Gauss_quadrature_1D_test(Neumann_fun,Gauss_nodes',Gauss_weights,Vertex_El,para.basis_type,alpha,0,0);
        end
    elseif boundary(1,k)==-3% Robin data
        El=boundary(2,k);% Element which as edge boundary(2,k).
        number_of_local_basis_test=size(Tb_test,1); %number of local basis function on element El
        number_of_local_basis_trial=size(Tb_trial,1);
        Vertex_El=Pb(:,Tb_test(:,El));%vertices of the element
        vertices=Pb(:,[boundary(3,k);boundary(4,k)]);%vertices of the edge
        edge_length=norm(vertices(:,1)-vertices(:,2),2);
        [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(vertices(:,1),vertices(:,2),para.order_Gauss);%already multiplied by length edge
        %Assemble additional contribution given by the integral on the
        %edge.
        for alpha=1:number_of_local_basis_trial %loop over trial functions on the element for trial FE space.
            for beta=1:number_of_local_basis_test %loop over trial functions on the element for test FE space.
                int_value=Gauss_quadrature_1D_trial_test(@(x,y) para.Robin,Gauss_nodes,Gauss_weights,Vertex_El,Vertex_El,para.basis_type,alpha,0,0,para.basis_type,beta,0,0);
                A(Tb_test(beta,El),Tb_trial(alpha,El))=A(Tb_test(beta,El),Tb_trial(alpha,El))+int_value;
            end
        end
        %Modify Right hand side
        for alpha=1:number_of_local_basis_test%Assemble right hand side
           b(Tb_test(alpha,El))=b(Tb_test(alpha,El))+ Gauss_quadrature_1D_test(Robin_fun,Gauss_nodes',Gauss_weights,Vertex_El,para.basis_type,alpha,0,0);
        end
        
    elseif boundary(1,k)==-1% Dirichlet Data
        index=boundary(3:end,k);% indexes of nodes which are on the edge boundary
       
        for j=1:length(index)
             v(index(j))=-1;
        end
    end
    
end
w=find(v~=-1);
b=b(w);
A=A(w,w);
end