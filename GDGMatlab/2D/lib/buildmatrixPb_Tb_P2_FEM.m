function [Pb,Tb]=buildmatrixPb_Tb_P2_FEM(P,T,Eb);
%buildmatrixPb_Tb_P2_FEM receives the mesh matrices P,T and the edge matrix
%Eb.
%It returns the matrices 
%Pb: contains location of FEM DOFs.
%Tb: is a 6xNumber_triangle matrix. Tb(1:3,k) coincides with T(1:3,k)while
%the Tb(4,k) is the middle point between the vertices T(1,k) and T(2,k).
% Tb(5,k) is the middle point between the vertices T(2,k) and T(3,k).
% Similarly for Tb(6,k)




% author: Tommaso Vanzan
%-------------------------------------------------------------------------
Pb=zeros(2,size(P,2)+size(Eb,2));%Pb has the size of P plus one middle point for each edge.
Pb(:,1:size(P,2))=P;
lastvertex=size(Pb,1);
Tb=zeros(6,size(T,2));
Tb(1:3,:)=T(1:3,:);
lastvertex=size(P,2);
for k=1:size(Eb,2)
T_left=Eb(3,k);%Triangle on the left edge
T_right=Eb(4,k);%Triangle on the right edge
Vertexini=Eb(1,k);%initial vertex segment
Vertexend=Eb(2,k);%final vertex segment
NewVertex_x=(P(1,Vertexini)+P(1,Vertexend))/2;
NewVertex_y=(P(2,Vertexini)+P(2,Vertexend))/2;
lastvertex=lastvertex+1;
Pb(:,lastvertex)=[NewVertex_x,NewVertex_y]'; %Insert new vertex in Pb
if(T_left~=-1)
    if(Vertexini==Tb(1,T_left) && Vertexend==Tb(2,T_left))
        Tb(4,T_left)=lastvertex;
    elseif(Vertexini==Tb(2,T_left) && Vertexend==Tb(3,T_left))
        Tb(5,T_left)=lastvertex;
     elseif(Vertexini==Tb(3,T_left) && Vertexend==Tb(1,T_left))
        Tb(6,T_left)=lastvertex;
    end
end
if(T_right~=-1)
    if(Vertexini==Tb(2,T_right) && Vertexend==Tb(1,T_right))
        Tb(4,T_right)=lastvertex;
    elseif(Vertexini==Tb(3,T_right) && Vertexend==Tb(2,T_right))
        Tb(5,T_right)=lastvertex;
     elseif(Vertexini==Tb(1,T_right) && Vertexend==Tb(3,T_right))
        Tb(6,T_right)=lastvertex;
    end
end
end