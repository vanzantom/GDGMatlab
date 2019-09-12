function result=FE_local_basis(x,vertices,basis_type,basis_index,der)
%FE_local_basis receives a point in x in which evaluete the basis function.
% basis_type stands for the kind of FE Space.
% basis_index is a variable containing the index of the local basis
% function. Ex: 1D linear I have two values.
% der tells me if I want the derivative or not.

%===================

%basis_type==%101 1D linear nodal basis functions
if basis_type==101
   xn=vertices(1);
   xnp1=vertices(2);
   h=xnp1-xn;
   if der==0
       if basis_index==1
              result=(xnp1-x)/h;
       elseif basis_index==2
              result=(x-xn)/h;
       else 
            fprintf('error basis_index');
       end 
   elseif der==1
       if basis_index==1
              result=-1/h;
       elseif basis_index==2
              result=1/h;
       else 
            fprintf('error basis_index');
       end
   elseif der>1
      result=0; 
   end
else
    fprintf('no such type stored');
end