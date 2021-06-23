function result=FE_local_basis_1D(x,v1,v2,basis_type,basis_index,der)
%FE_local_basis_1D receives a point in x in which evaluete the basis function.
% basis_type stands for the kind of FE Space.
% basis_index is a variable containing the index of the local basis
% function. Ex: 1D linear I have two values.
% der tells me if I want the derivative or not.

%===================

%basis_type==%2010 1D linear nodal basis functions on a edge in 2D
if basis_type==2010
   if der==0
       if basis_index==1
              result=1-norm(x-v1)/norm(v2-v1);
       elseif basis_index==2
              result=norm(x-v1)/norm(v2-v1);
       else 
            fprintf('error basis_index');
       end 
   elseif der>0
      result=0;
      fprintf('no such type stored');
   end
else
   
end