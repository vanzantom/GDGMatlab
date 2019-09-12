function [P,T,E,Pb,Tb,Eb,h]=generate_mesh_2D(h,basis_type)
% The function 'generate_mesh_2D' generates the mesh on a geometry specified by a geometric function with
% maximum mesh size 'h'. It returns the vertex matrix P, the element matrix
% T associated to the mesh and the vertex and Element matrix corresponding
% to the basis_type.
% Structure T(:,1)=[#V_1;#V_2;#V_3;Nsubdomain].
%           P(:,1)=[X_1;Y_1]
%In E, the first and second rows contain indices of the starting and ending point, 
%the third and fourth rows contain the starting and ending parameter values, 
%the fifth row contains the edge segment number, and the sixth and seventh row contain 
%the left- and right-hand side subdomain numbers.
%===================================================
% Example:
%201: 2D linear nodal FE on triangle
%===== Generation Mesh
[P,E,T]=initmesh(@square_twosub,'Hmax',h);%general triangular mesh
%[P,E,T] = poimesh(@square_twosub,1/h,1/h);%regular triangular mesh
%==== Data structure for P1 FE
if basis_type==201
    Pb=P;
    Tb=T(1:end-1,:);
    Eb=E;%not needed.
%==== Data structure for DG P1 FE
% Output: Eb the first and second rows contain indices of the starting and
% ending vertices
% the third and fourth rows contains triangle on the left and on the right
% according to anticlockwise orientation
% Output: T has 7 rows. 3 Vertices 1 Subdomain 3 Edges
elseif basis_type==2010
    Tb=T(1:end-1,:);
    T=[T;zeros(3,size(T,2))]; 
    j=0;
    for k=1:size(T,2) % loop over the element
          for i=1:3 % I assigned the vertices of the DOF to the vertices of the elements.
              j=j+1;
              Pb(:,j)=P(:,T(i,k));
          end
          Tb((1:3),k)=[j-2;j-1;j]; % Now each Triangule has its own DOF.
    end
    %% ============= BUILD MATRIX E_B
max_num_edges= 3 * size(P,2)- 6;%estimate to preallocated a sufficiently large matrix Eb
num_edges=1;
Eb=-ones(4,max_num_edges); %Initial all right triangles to the edge as "boundary". Then if it is not the case we change the value, inserint T_right
vertex1=T(1,1);%First triangle, I insert all edges.
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
Eb(3,num_edges)=1;
T(5,1)=1;T(6,1)=2;T(7,1)=3;% I insert the edges label in the matrix T.
for k=2:size(T,2)%new triangle => new 3 edges.
    vertex1=T(1,k);%save vertices.
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
    
    
%===== Compute diameter mesh. h=minimum length of an edge.
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