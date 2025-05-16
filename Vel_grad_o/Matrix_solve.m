% x_mat = mldivide(MatA_prsd,MatB_prsd);
% 

x_mat_raw = mldivide(MatA_prsd,MatB_prsd);
x_mat=zeros(3*N_total+Nx*Ny,1);
for n = 1:3*N_total+Nx*Ny
    if n < 3*N_total
        x_mat(n) = x_mat_raw(n);
    elseif n > 3*N_total
        x_mat(n) = x_mat_raw(n-1);
    end
end

R_mat = zeros(N_total,3);

for n = 1:3*N_total
    remainder = rem(n,3);
    if (remainder == 0)
        R_mat(n/3 ,3) = x_mat(n);
    elseif (remainder == 1)
        R_mat((n+2)/3,1) = x_mat(n);
    elseif (remainder == 2)
        R_mat((n+1)/3,2) = x_mat(n);
    end

end


for n = 1:Nx*Ny
    Velocity_boundaries(n+Nx*Ny+Nx*Nz,3) = x_mat(3*N_total + n);
end