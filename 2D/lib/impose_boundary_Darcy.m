function [A,ARobin,ADir,b,M_GG,e]=impose_boundary_Darcy(ANeumann,b,boundary,Pb,Tb_trial,Tb_test,para,Dirichlet_fun,Neumann_fun)
%-------------------------------------------------------------------------
% Impose_boundary_Darcy receives 
% A: stiffness matrix
% b: right hand side vector
% boundary: boundaryedge matrix
% Pb: matrix with position degrees of freedom.
% Tb_trial: matrix with elements information for trial FE space
% Tb_test: matrix with elements information for test FE space
% Dirichlet_fun: Dirichlet boundary data
% Neumann_fun: Dirichlet boundary data
%treat_boundary returns
% A: modified stifness matrix with Dirichlet BC in -1 and Neumann
% for -2 and -1
% ADir: modified stiffness matrix where the Dirichlet BC are strongly
% imposed on -1 and -3 (needed for DtN)
% ARobin: modified stiffness matrix where additional Robin BC are imposed
% on -3
% b: modified right hand side vector to include Dirichlet and Neumann boundary conditions
% M_GG: local mass matrix which contains only the mass terms on the interface.
% vector e: containing index vertices on the edge
% author: Tommaso Vanzan
%-------------------------------------------------------------------------
%========================================================================
nbn=size(boundary,2);
tol=10^-8;
Mass=sparse(size(ANeumann,1),size(ANeumann,2));
%=======================
% Find indexes nodes on the interface!
%======================
l=1;
for k=1:size(Pb,2)
    if(Pb(1,k)>tol && Pb(1,k)<1-tol && abs(Pb(2,k))<tol)
         e(l)=k;
         l=l+1;
    end
end
[~,i]=sort(Pb(1,e));%Reorder vector e such that interface nodes are ordered!
esup=e;
for l=1:length(i)
    e(l)=esup(i(l));
end



%=============
% Assemble Neumann BC (if needed) in b and MAss matrix on the interface
%=============


Mass=sparse(size(ANeumann,1),size(ANeumann,2));;
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
                int_value=Gauss_quadrature_1D_trial_test(@(x,y) 1,Gauss_nodes,Gauss_weights,Vertex_El,Vertex_El,para.basis_type,alpha,0,0,para.basis_type,beta,0,0);
                Mass(Tb_test(beta,El),Tb_trial(alpha,El))=Mass(Tb_test(beta,El),Tb_trial(alpha,El))+int_value;
            end
        end
    end
end

M_GG= Mass(e,e);% Restricted mass matrix
A=ANeumann;%Impose Dirichlet BC on the Neumann Matrix
for k=1:nbn
    if boundary(1,k)==-1% Dirichlet Data
        index=boundary(3:end,k);% indexes of nodes which are on the edge boundary
        for j=1:length(index)
            A(index(j),:)=0;
            A(index(j),index(j))=1;
            b(index(j))=feval(Dirichlet_fun,Pb(1,index(j)),Pb(2,index(j)));
        end
    end
end

ARobin=A;% Add Robin term!
for i=1:length(e)
    ARobin(e(i),:)=ARobin(e(i),:)+1/para.Robin*Mass(e(i),:);
end


ADir=ANeumann;
for k=1:nbn
    if(boundary(1,k)==-1  || boundary(1,k)==-3)
       index=boundary(3:end,k);% indexes of nodes which are on the edge boundary
       for j=1:length(index)
            ADir(index(j),:)=0;
            ADir(index(j),index(j))=1;
       end
    end
end


end

