%-------------------------------------------------------------------------
% treat_boundary receives 
% A: stiffness matrix
% b: right hand side vector
% boundaryu/boundaryv: boundary data structure for the two velocities
% P: data structure for vertex mesh
% Pb: matrix with position degrees of freedom.
% Tb_test: data structure DOF velocity
% Dirichlet_fun: Dirichlet boundary data
%treat_boundary returns
% ADir: modified stiffness matrix where the Dirichlet BC are strongly
% imposed on every boundary edge -1 or -3 (needed to build DtN)
% A: modified stifness matrix with Dirichlet BC in -1 and Neumann
% for -2 and -3
% ARobin: includes also the mass term on the interface labeled by -3.
% b: modified right hand side vector to include Dirichlet and Neumann boundary
% conditions
% vector e: containing index vertices on the edge
% author: Tommaso Vanzan
%-------------------------------------------------------------------------

function [A,ARobin,ADir,b,e,M_GG]=impose_boundary_Stokes(ANeumann,b,boundaryu,boundaryv,P,Pb,Tb_test,Tb_trial,data,para)

number_of_local_basis_trial=size(Tb_trial,1);
number_of_local_basis_test=size(Tb_test,1);
%=======================
% Find indexes nodes on the interface!
%======================
ndofve=size(Pb,2);
l=1;
tol=10^(-8);
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
e=e+ndofve;%increment indexes by ndofve to find DOF on interface for velocity v
%========================
% Create Mass matrix
%========================
nbn=size(boundaryu,2);
A=ANeumann;
Mass=sparse(size(A,1),size(A,2));
for k=1:nbn
    if boundaryv(1,k)==-2
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
            El=boundaryv(2,k);% Element which as edge boundary(2,k).
            number_of_local_basis_test=size(Tb_test,1); %number of local basis function on element El
            Vertex_El=P(:,Tb_test(1:3,El));%vertices of the element
            vertices=P(:,[boundaryv(3,k);boundaryv(4,k)]);%vertices of the edge
            edge_length=norm(vertices(:,1)-vertices(:,2),2);
            %normal=find_normal(boundaryv(3,k),boundaryv(4,k),Tb_test(1:3,El),P);
            [Gauss_nodes,Gauss_weights]=generate_Gauss_1D(vertices(:,1),vertices(:,2),para.order_Gauss);%already multiplied by length edge
                for alpha=1:number_of_local_basis_trial %loop over trial functions on the element for trial FE space.
                    for beta=1:number_of_local_basis_test %loop over trial functions on the element for test FE space.
                        int_value=Gauss_quadrature_1D_trial_test(@(x,y) 1,Gauss_nodes,Gauss_weights,Vertex_El,Vertex_El,para.basis_type_ve,alpha,0,0,para.basis_type_ve,beta,0,0);
                        Mass(Tb_test(beta,El)+ndofve,Tb_trial(alpha,El)+ndofve)=Mass(Tb_test(beta,El)+ndofve,Tb_trial(alpha,El)+ndofve)+int_value;
                    end
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
ARobin=A;
for i=1:length(e)
    ARobin(e(i),:)=ARobin(e(i),:)+para.Robin*Mass(e(i),:);
end
M_GG= Mass(e,e);% Restricted mass matrix

%========================
% Create a matrix with Dirichlet everywhere
%========================
% ADir=ANeumann;
% for k=1:nbn
%     index=boundaryu(3:end,k);% indexes of nodes which are on the edge boundary
%         for j=1:length(index)
%             ADir(index(j),:)=0;
%             ADir(index(j),index(j))=1;
%         end
%     
%        index=boundaryv(3:end,k);% indexes of nodes which are on the edge boundary
%        for j=1:length(index)
%             ADir(index(j)+ndofve,:)=0;
%             ADir(index(j)+ndofve,index(j)+ndofve)=1;
%        end
% end

%========================
% Create a matrix with Dirichlet where it is supposed to be plus where
% Robin data should stay
%========================
ADir=ANeumann;
for k=1:nbn
    if(boundaryu(1,k)==-1)
    index=boundaryu(3:end,k);% indexes of nodes which are on the edge boundary
        for j=1:length(index)
            ADir(index(j),:)=0;
            ADir(index(j),index(j))=1;
        end
    end
    if(boundaryv(1,k)==-1  || boundaryv(1,k)==-3)
       index=boundaryv(3:end,k);% indexes of nodes which are on the edge boundary
       for j=1:length(index)
            ADir(index(j)+ndofve,:)=0;
            ADir(index(j)+ndofve,index(j)+ndofve)=1;
       end
    end
end


end



    