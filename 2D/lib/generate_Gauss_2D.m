function [Gauss_nodes,Gauss_weights]=generate_Gauss_2D(vertices,type)
% The function 'generate_Gauss_2D' delivers the Gaussian nodes on the
% current element defined by the vector vertices.


%=== Linear map 
A=[(vertices(1,2)-vertices(1,1)),(vertices(1,3)-vertices(1,1));(vertices(2,2)-vertices(2,1)),(vertices(2,3)-vertices(2,1))];

if type >=1
    %=== Order of Gauss quadrature and nodes on the reference element.
    w_2D=[];
    node_2D=[];
    % nodes on the 1d interval [-1,1]
    n=type;
    [x,w]= gauleg(-1,1,n);
    for i=1:n
        for j=1:n
            node_2D=[node_2D; (1+x(i))./2 , (1-x(i)).*(1+x(j))./4];
            w_2D = [w_2D, (1-x(i)).*w(i).*w(j)./8];
        end
    end
    v=[vertices(1,1); vertices(2,1)];
    v=repmat(v,1,size(node_2D,1));
    Gauss_nodes=node_2D;
    Gauss_nodes=A*Gauss_nodes'+v;
    Gauss_nodes=Gauss_nodes';
    Gauss_weights=w_2D;
elseif type==-1
    Gauss_nodes=vertices';
    Gauss_weights=[1/3,1/3,1/3];
elseif type==-2   
    %Lumped for P2b 
    Gauss_nodes=[0,0;1,0;0,1;0.5,0;0,0.5;0.5,0.5;1/3,1/3];
    v=[vertices(1,1); vertices(2,1)];
    v=[v,v,v,v,v,v,v];
    Gauss_nodes=A*Gauss_nodes'+v;
    Gauss_nodes=Gauss_nodes';
    Gauss_weights=[1/20,1/20,1/20,2/15,2/15,2/15,9/20];
end
    
    
function [x,w]= gauleg(a,b,n)

m=(n+1)/2;
xm=0.5*(b+a);
xl=0.5*(b-a);
xx=[];

for i=1:m
    z=cos(pi*(i-0.25)/(n+0.5));
    while 1
        p1=1.0;
        p2=0.0;
        for j=1:n
            p3=p2;
            p2=p1;
            p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
        end
        pp=n*(z*p1-p2)/(z*z-1.0);
        z1=z;
        z=z1-p1/pp;
        if (abs(z-z1)<eps), break, end
    end
    xx(i)=xm-xl*z;
    xx(n+1-i)=xm+xl*z;
    ww(i)=2.0*xl/((1.0-z*z)*pp*pp);
    ww(n+1-i)=ww(i);
end

x=xx;
w=ww;
end

end
%Newton Cotes
% if type==1
% Gauss_nodes(1,1)=1/3; Gauss_nodes(1,2)=1/3;
% Gauss_weights(1)=1;
% elseif type==2
% Gauss_nodes(1,1)=0;Gauss_nodes(1,2)=1/2;
% Gauss_nodes(2,1)=1/2; Gauss_nodes(2,2)=0;
% Gauss_nodes(3,1)=1/2; Gauss_nodes(3,2)=1/2;
% Gauss_weights(1)=1/3;
% Gauss_weights(2)=1/3;
% Gauss_weights(3)=1/3;
% end
% %=== Map node on the current element.
% v=[vertices(1,1); vertices(2,1)];
% v=[v,v,v];
% Gauss_nodes=A*Gauss_nodes'+v;
% Gauss_nodes=Gauss_nodes';
% Gauss_weights=Gauss_weights/2;
% end