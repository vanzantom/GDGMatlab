%-------------------------------------------------------------------------
% generate_boundaryedges receives matrix Eb, the mesh matrix E and the
% structures data and para.
% generate_boundaryedges returns: 
% boundary (matrix
% boundary(1,j): starting vertex boundary edge j
% boundary(2,j): ending  vertex boundary edge  j
% boundary(3,j): element to which boundary element j belongs
% boundary(4,j): label of the boundary(Dirichlet -1,Neumann -2,Robin -3)
% boundary(5,j): possible middle point for P2 FEM.


% author: Tommaso Vanzan
%-------------------------------------------------------------------------

function boundary=generate_boundaryedges(Tb,Eb,E,label,basis_type)
v=find(Eb(4,:)==-1); %indexes of edges which have -1 as Right element.
boundary=zeros(4,length(v));
flag4=0;
boundary(3,:)=Eb(1,v);
boundary(4,:)=Eb(2,v);
boundary(2,:)=Eb(3,v);
for i=1:length(v)
    j=1;
    while j<=size(E,2) && flag4==0
        if(sum(Eb(1:2,v(i))==E(1:2,j))==2 || sum(Eb([2,1],v(i))==E(1:2,j))==2) %check to which edge in E, Eb(,v(i)) corresponds!
            index=E(5,j);
            flag4=1;
        end
        j=j+1;
    end
    boundary(1,i)=label(index);%assign label given to edge index in E!
    flag4=0;
end

if(basis_type==202)%P2 FE. Insert middle point of the edge.
    boundary=[boundary;zeros(1,size(boundary,2))];
    for i=1:size(boundary,2)
        El=boundary(2,i);
        if(sum(boundary(3:4,i)==Tb(1:2,El))  )
            boundary(5,i)=Tb(4,El);
        elseif(sum(boundary(3:4,i)==Tb(2:3,El))  )
            boundary(5,i)=Tb(5,El);
        elseif(sum(boundary(3:4,i)==Tb([3,1],El))  )
            boundary(5,i)=Tb(6,El);
        end
    end
end
        
end
