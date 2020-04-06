%-------------------------------------------------------------------------
% The function OSM performs one iteration of optimized Schwarz method for
% HDG discretization, written at the substructured level.
% OSM receives:
% lambdas: list of Nsubdomains vectors, lambdas{i} is the substructured
%           related to subdomain i
% AII, AGamma, AIGAMMa: list of matrices
% B: list of matrix in case I use modified OSM using a relaxation parameter
% R: restriction operators on each subdomain structure.
% Nsubdomains: number of subdomains.
% gamma: relaxation parameter.

% author: Tommaso Vanzan
%-------------------------------------------------------------------------


function [lambdasnew]=OSM(lambdas,Aii,AGamma,AIGamma,B,R,b,Nsubdomains,gamma)

lambdasnew=cell(Nsubdomains,1);
for i=1:Nsubdomains
        f{i}=R{i}*b;
        if(gamma==0) %if gamma==0, standard OSM for HDG
        for k=1:Nsubdomains
            if(k~=i)
                f{i}=f{i}+R{i}*AIGamma{k}'*(Aii{k}\(AIGamma{k}*R{k}'*lambdas{k}));
            end
        end
        lambdasnew{i}=(AGamma{i}-B{i})\f{i};
        f{i}=0*f{i};
        end
        if(gamma~=0) %if gamma==0, I perform modified OSM for HDG discretizations
         for k=1:Nsubdomains
            if(k~=i)
                f{i}=f{i}+R{i}*AIGamma{k}'*(Aii{k}\(AIGamma{k}*R{k}'*lambdas{k})) -(1-gamma)*AGamma{k}*lambdas{k};
            end
         end
            lambdasnew{i}=(gamma*AGamma{i}-B{i})\f{i};
            f{i}=0*f{i};   
         end
    end
end
