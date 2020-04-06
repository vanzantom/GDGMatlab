function [P,T,E,Pb,Tb,Eb,h,Nsub,ngamma,indexGamma]=generate_mesh_2D_HDG_backup(h,basis_type)
% The function 'generate_mesh_2D_HDG'  generates the mesh on a geometry specified by a geometric function with
% maximum mesh size 'h'. It returns the vertex matrix P, the element matrix
% T associated to the mesh and the vertex and Element matrix corresponding
% to the basis_type.
% Structure T(:,1)=[#V_1;#V_2;#V_3;Nsubdomain].
%           P(:,1)=[X_1;Y_1]
%In E, the first and second rows contain indices of the starting and ending point, 
%the third and fourth rows contain the starting and ending parameter values, 
%the fifth row contains the edge segment number, and the sixth and seventh row contain 
%the left- and right-hand side subdomain numbers.
%It returns also the number of DOF in first and second subdomain
%nsub1,nsub2 and the number of DOF on the interface gamma.
% Example:
%===============================
%201: 2D linear nodal FE on triangle
%===== Generation Mesh
[P,E,T]=initmesh(@square_foursub,'Hmax',h);%general triangular mesh
%[P,E,T]=initmesh(@square_twosub,'Hmax',h);%general triangular mesh
%[P,E,T] = poimesh(@square_twosub,1/h,1/h);%regular triangular mesh. NOT
%compatible with HDG.
%==== Data structure for P1 FE
if basis_type==201
    Pb=P;
    Tb=T(1:end-1,:);
    Eb=E;%not needed.
%==== Data structure for HDG P1 FE
% Output: Eb the first and second rows contain indices of first and second
% vertex
% the third and fourth rows contains triangle on the left and on the right
% according to anticlockwise orientation. 5th rows equal to 1 if the edge
% lies on the interface between subdomains. 6th and 7th rows contains the
% number of the associated degree of freedoms on the edge on the interface
% Output: T has 7 rows. 3 index Vertices, 1 Subdomain label, 3 Edges
% indexes
% parameters nsub1 and nsub2 describe last DOF index of subdomain 1 and 2.
elseif basis_type==2010
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
j=Nsub(end);
indeces=ones(Nsubdomains);%index vectors. One counting variable for each subdomain
indexGamma=cell(Nsubdomains,1);% Each cell contains the indeces of the DOF on GAMMA belonging to each subdomain
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
            indexGamma{T_1_label}(indeces(T_1_label))=j-1;%IndexGamma{1} is a vector containing the DOFs on Gamma which lie on Gamma_1
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