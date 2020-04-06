function [lambdasnew]=OSM_damped(lambdas,Aii,AGamma,AIGamma,B,R,b,Nsubdomains,w)

lambdasnew=cell(Nsubdomains,1);
for i=1:Nsubdomains
        f{i}=R{i}*b;
        for k=1:Nsubdomains
            if(k~=i)
                f{i}=f{i}+R{i}*AIGamma{k}'*(Aii{k}\(AIGamma{k}*R{k}'*lambdas{k}));
            end
        end
        lambdasnew{i}=w*( ((1-w)/w)*lambdas{i} +(AGamma{i}-B{i})\f{i});
        f{i}=0*f{i};
        
end
end
