function normal=find_normal(vertex1,vertex2,vertex_left_index,P)
% The function find_normal receives two indeces of two vertices and a
% vector of three vertices indeces. It finds the third missing index
% vertex3. Then it computes the tangential vector tau to the edge, the
% normal one and computing the scalar product btw the tau and normal
% decides the sign of the normal vector.
normal=zeros(2,1);
if(vertex1==vertex_left_index(1) && vertex2==vertex_left_index(2) || vertex2==vertex_left_index(1) && vertex1==vertex_left_index(2))
         vertex3=vertex_left_index(3);
    elseif(vertex1==vertex_left_index(2) && vertex2==vertex_left_index(3) || vertex2==vertex_left_index(2) && vertex1==vertex_left_index(3)) 
         vertex3=vertex_left_index(1);
    elseif(vertex1==vertex_left_index(3) && vertex2==vertex_left_index(1) || vertex2==vertex_left_index(3) && vertex1==vertex_left_index(1)) 
         vertex3=vertex_left_index(2);
end
    tau=P(:,vertex2)-P(:,vertex1); %compute tangential vector
    tau=tau/norm(tau);%normalize tangential vector
    normal(1)=-tau(2);
    normal(2)=tau(1);
    vintern=P(:,vertex3)-P(:,vertex1);
    vintern=vintern/norm(vintern);
    if vintern'*normal >0
        normal=-normal;
    end
end