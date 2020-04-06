function [Eb,Pb,indexGamma,j]=order_vertices_Gamma_twosubdomains(Eb,T,P,Pb,indexGamma,indeces,num_edges,j)
count=1;
V=[];
    for k=1:num_edges
        T_1= Eb(3,k);
        T_2=Eb(4,k);
        if(T_2~=-1)
            T_1_label=T(4,T_1);
            T_2_label=T(4,T_2);
            if(T_1_label~=T_2_label)
                V(1:2,count)=P(:,Eb(1,k));%read edge's vertices.
                V(3:4,count)=P(:,Eb(2,k));
                if(V(2,count)>V(4,count)) %I order the edge from the bottom to the top.
                    temp=V(1:2,count);
                    V(1:2,count)=V(3:4,count);
                    V(3:4,count)=temp;
                    temp=Eb(1,k);
                    Eb(1,k)=Eb(2,k);
                    Eb(2,k)=temp;
                    temp=Eb(3,k);
                    Eb(3,k)=Eb(4,k);
                    Eb(4,k)=temp;
                end
                V(5,count)=k;
                V(6,count)=T_1_label;
                V(7,count)=T_2_label;
                for i=1:size(V,2)-1  % force order between new vertices and the one already inserted.
                    if(V(2,count)<=V(2,i) && V(4,count)<=V(4,i))
                        temp=V(:,i);
                        V(:,i)=V(:,count);
                        V(:,count)=temp;
                    end
                end
                count=count+1;
            end
        end
    end %At the end of this loop, V contains all the information I need while still having the vertices ordered.
    for i=1:size(V,2)
        Eb(5,V(5,i))=1;
        j=j+1;
        Pb(:,j)=V(1:2,i);
        j=j+1;
        Pb(:,j)=V(3:4,i);
        indexGamma{V(6,i)}(indeces(V(6,i)))=j-1;%IndexGamma{k} is a vector containing the DOFs on Gamma which lie on Gamma_k(Gamma_k is the interface of subdomain k)
        indexGamma{V(6,i)}(indeces(V(6,i))+1)=j;
        indexGamma{V(7,i)}(indeces(V(7,i)))=j-1;
        indexGamma{V(7,i)}(indeces(V(7,i))+1)=j;
        indeces(V(6,i))=indeces(V(6,i))+2;
        indeces(V(7,i))=indeces(V(7,i))+2;
        Eb(6,V(5,i))=j-1;
        Eb(7,V(5,i))=j;
    end
end