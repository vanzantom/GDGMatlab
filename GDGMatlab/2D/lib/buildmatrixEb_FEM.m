function [Eb]=buildmatrixEb_FEM(P,T);
%buildmatrixE_FEM creates a matrix E of edges.
%It receives
% P: matrix of mesh nodes
% T: matrix of triangles
% It returns
%matrix E whose rows are 
%E(:,1) starting vertex
%E(:,2) ending vertex
%E(:,3) triangle on the left going from E(:,1) to E(:,2)
%E(:,4) triangle on the right going from E(:,1) to E(:,2)

% author: Tommaso Vanzan
%-------------------------------------------------------------------------
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
        flag_1=1;
    elseif ( (vertex2==Eb(2,j) && vertex3==Eb(1,j)) || (vertex2==Eb(1,j) && vertex3==Eb(2,j)) ) %check if second edge is already inside.
        Eb(4,j)=k;
        flag_2=1;
    elseif ( (vertex3==Eb(2,j) && vertex1==Eb(1,j)) || (vertex3==Eb(1,j) && vertex1==Eb(2,j))) %check if third edge is already inside.
        Eb(4,j)=k;
        flag_3=1;
    end
end
if(flag_1==0) %If there are not inside, then I add them.
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex1;
    Eb(2,num_edges)=vertex2;
    Eb(3,num_edges)=k;
end
if(flag_2==0)
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex2;
    Eb(2,num_edges)=vertex3;
    Eb(3,num_edges)=k;
end
if(flag_3==0)
    num_edges=num_edges+1;
    Eb(1,num_edges)=vertex3;
    Eb(2,num_edges)=vertex1;
    Eb(3,num_edges)=k;
end

end
Eb=Eb(:,(1:num_edges)); 


end