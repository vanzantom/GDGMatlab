function result=FE_local_basis_2D(x,y,vertices,basis_type,basis_index,der_x,der_y)



J=abs((vertices(1,2)-vertices(1,1))*(vertices(2,3)-vertices(2,1)) -(vertices(1,3)-vertices(1,1))*(vertices(2,2)-vertices(2,1)));

x_hat=((vertices(2,3)-vertices(2,1))*(x-vertices(1,1))-(vertices(1,3)-vertices(1,1))*(y-vertices(2,1)))/J;
y_hat=((vertices(2,1)-vertices(2,2))*(x-vertices(1,1))-(vertices(1,1)-vertices(1,2))*(y-vertices(2,1)))/J;


if der_x==0 && der_y==0
    result=FE_reference_basis_2D(x_hat,y_hat,basis_type,basis_index,der_x,der_y);
elseif der_x==1 && der_y==0
    result= ((vertices(2,3)-vertices(2,1)))/J*FE_reference_basis_2D(x_hat,y_hat,basis_type,basis_index,1,0) +1/J*(-vertices(2,2)+vertices(2,1))*FE_reference_basis_2D(x_hat,y_hat,basis_type,basis_index,0,1);  
elseif der_x==0 && der_y==1
    result= ((-vertices(1,3)+vertices(1,1)))/J*FE_reference_basis_2D(x_hat,y_hat,basis_type,basis_index,1,0) + 1/J*(vertices(1,2)-vertices(1,1))*FE_reference_basis_2D(x_hat,y_hat,basis_type,basis_index,0,1);
end


end