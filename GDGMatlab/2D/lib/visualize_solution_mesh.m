function []=visualize_solution_mesh(T,Pb,Tb,u,flag)
% The function visualize_mesh receives the mesh matrices Pb,Tb and T and the
% number of interior degree of freedoms n and a flag
% if flag==0 it plots the solution
% if flag==1, it shows the mesh elements. If flag==2, it shows the subdomain
% decomposition
n=length(u);
Tb=[Tb;T(4,:)];
Pb=Pb(:,1:n);
if(flag==0)
    pdeplot(Pb,[],Tb,'xydata',u,'xystyle','flat',...
           'zdata',u,'zstyle','continuous','colorbar','on');
elseif(flag==1)
u=zeros(n,1);
for i=1:size(Tb,2)
   u(Tb((1:3),i))=i;
end
elseif(flag==2)
u=zeros(n,1);

for i=1:size(Tb,2)
    if(Tb(4,i)==1)
    u(Tb((1:3),i))=1;
    end
    if(Tb(4,i)==2)
    u(Tb((1:3),i))=2;
 end
end
end

pdeplot(Pb,[],Tb,'xydata',u,'xystyle','flat',...
           'zdata',u,'zstyle','continuous','colorbar','on');
end