%-------------------------------------------------------------------------
% generate_mesh_2D_HDG generates the mesh on a geometry specified by the geometric function geo.fun with
% maximum mesh size geo.h.
% generate_mesh_2D_HDG receives
% geo: structure containing a function handle describing the geometry and
% the meshsize geo.h and a flag for a regular/non regular mesh.
% basis_type: type of DG space (linear,quadratic etc..).
% generate_mesh_2D_HDG returns:
% P: vertex matrix containing position of vertices of the mesh. Eg: P(:,1)=[X_1;Y_1]
% T: element matrix containing indeces vertices of the elements plus subdomain to which element k belongs. 
%     Eg: T(:,1)=[#V_1;#V_2;#V_3;Nsubdomain]. If using DG, in T we also
%     save the indexes of the edges of the triangles.
% E: edge matrix. The first and second rows contain indices of the starting and ending point of the edge, 
%   the third and fourth rows contain the starting and ending parameter values, 
%   the fifth row contains the edge segment number, and the sixth and seventh row contain 
%   the indeces of subdomains on the left and on the right side.
% Pb: vertex matrix of the degrees of freedom (If using P1, Pb=P).
% Tb: element matrix with indeces degrees of freedom on each element ((If using P1, Tb=T)
% Eb: edge matrix(not needed for FEM, only for DG!). In Eb the first and second rows contain indices of the starting and
% ending vertices
% the third and fourth rows contains triangle on the left and on the right
% according to anticlockwise orientation
% h: a scalar variable which is equal to the minimum length of an edge.
% Nsub:   vector. Nsub(j) contains last index DoF in subdomain j.
% ngamma:  total number of DOF on Gamma
% indexGamma: vector containing indexes DOFs on Gamma. (If using two
%               subdomains, the indexes in indexGamma are ordered on the
%               interfaces increasingly in y

% author: Tommaso Vanzan
%-------------------------------------------------------------------------


function [P,T,E,Pb,Tb,Eb,h,Nsub,ngamma,indexGamma]=generate_mesh_2D_HDG(geo,basis_type)
    

[P,E,T]=initmesh(geo.fun,'Hmax',geo.h);%general triangular mesh

if basis_type==2010
    Tb=T(1:end-1,:);% Just line to copy rows of T.
    T=[T;zeros(3,size(T,2))];% 
    j=0;
    %==== I run over the subdomains and triangles to assign DOF in order to
    %each subdomain!
    Nsubdomains=max(T(4,:)); %get Numbers of subdomains
    for jj=1:Nsubdomains
        for k=1:size(T,2) % loop over the elements
            if(T(4,k)==jj)
              for i=1:3 % I assigned the vertices of the DOF to the vertices of the elements.
                  j=j+1;
                  Pb(:,j)=P(:,T(i,k));
              end
              Tb((1:3),k)=[j-2;j-1;j]; % Now each Triangule has its own DOF.
            end
        end
        Nsub(jj)=j;%last index is the number of DOFs in Nsub
    end
    %% ============= BUILD MATRIX E_B
max_num_edges= 3 * size(P,2)- 6;
num_edges=1;
Eb=-ones(4,max_num_edges); %We put edges from first triangule and relative edges.
vertex1=T(1,1);
vertex2=T(2,1);
vertex3=T(3,1);
Eb(1,num_edges)=vertex1;
Eb(2,num_edges)=vertex2;
Eb(3,num_edges)=1;
num_edges=num_edges+1;
Eb(1,num_edges)=vertex2;
Eb(2,num_edges)=vertex3;
Eb(3,num_edges)=1;
num_edges=num_edges+1;
Eb(1,num_edges)=vertex3;
Eb(2,num_edges)=vertex1;
Eb(3,num_edges)=1; % I inserted first three edges
T(5,1)=1;T(6,1)=2;T(7,1)=3;
for k=2:size(T,2)%new triangle => new 3 edges.
    vertex1=T(1,k);
    vertex2=T(2,k);
    vertex3=T(3,k);
    m=num_edges;
    flag_1=0;
    flag_2=0;
    flag_3=0;
for j=1:m
    if( (vertex1==Eb(2,j) && vertex2==Eb(1,j)) || (vertex1==Eb(1,j) && vertex2==Eb(2,j)) ) %check if first edge is already inside.
        Eb(4,j)=k;% I put the neibouring triangle to Eb equal to k
        T(5,k)=j; % I put the edge in the triangle matrix.
        flag_1=1;
    elseif ( (vertex2==Eb(2,j) && vertex3==Eb(1,j)) || (vertex2==Eb(1,j) && vertex3==Eb(2,j)) ) %check if second edge is already inside.
        Eb(4,j)=k;
        T(6,k)=j;
        flag_2=1;
    elseif ( (vertex3==Eb(2,j) && vertex1==Eb(1,j)) || (vertex3==Eb(1,j) && vertex1==Eb(2,j))) %check if third edge is already inside.
        Eb(4,j)=k;
        T(7,k)=j;
        flag_3=1;
    end
end
if(flag_1==0) %If there are not inside, then I add them.
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex1;
    Eb(2,num_edges)=vertex2;
    Eb(3,num_edges)=k;
    T(5,k)=num_edges;
end
if(flag_2==0)
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex2;
    Eb(2,num_edges)=vertex3;
    Eb(3,num_edges)=k;
    T(6,k)=num_edges;
end
if(flag_3==0)
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex3;
    Eb(2,num_edges)=vertex1;
    Eb(3,num_edges)=k;
    T(7,k)=num_edges;
end

end
Eb=Eb(:,(1:num_edges));   
Eb=[Eb;zeros(3,size(Eb,2))];% interface label
%========================================
% Identify which edges are on the interface and add DOF in Pb.
%========================================
j=Nsub(end);
indeces=ones(Nsubdomains);%index vectors. One counting variable for each subdomain
indexGamma=cell(Nsubdomains,1);% Each cell contains the indeces of the DOF on GAMMA belonging to each subdomain
if(strcmp(geo.fun,'square_twosub')==0) %Unless I am dealing with two subdomains, I do not care about the order of the points on the interface.
    for k=1:num_edges
        T_1= Eb(3,k);
        T_2=Eb(4,k);
        if(T_2~=-1)
            T_1_label=T(4,T_1);
            T_2_label=T(4,T_2);
            if(T_1_label~=T_2_label);
                    Eb(5,k)=1;
                    j=j+1;
                    Pb(:,j)=P(:,Eb(1,k));
                    j=j+1;
                    Pb(:,j)=P(:,Eb(2,k));
                    indexGamma{T_1_label}(indeces(T_1_label))=j-1;%IndexGamma{k} is a vector containing the DOFs on Gamma which lie on Gamma_k(Gamma_k is the interface of subdomain k)
                    indexGamma{T_1_label}(indeces(T_1_label)+1)=j;
                    indexGamma{T_2_label}(indeces(T_2_label))=j-1;
                    indexGamma{T_2_label}(indeces(T_2_label)+1)=j;
                    indeces(T_1_label)=indeces(T_1_label)+2;
                    indeces(T_2_label)=indeces(T_2_label)+2;
                    Eb(6,k)=j-1;
                    Eb(7,k)=j;
            end
        end
    end
end
if(strcmp(geo.fun,'square_twosub')==1)
   [Eb,Pb,indexGamma,j]=order_vertices_Gamma_twosubdomains(Eb,T,P,Pb,indexGamma,indeces,num_edges,j);
end
    
    
    
    
ngamma=j-Nsub(end);

%% Compute diameter mesh
h=0;
for i=1:num_edges
vertex1=Eb(1,i);
vertex2=Eb(2,i);
vertices1=P(:,vertex1);
vertices2=P(:,vertex2);
edge_length=norm(vertices1-vertices2,2);
if (edge_length>h)
   h=edge_length;
end
end
end
end