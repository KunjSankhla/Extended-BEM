MatA_prsd = MatA;
MatB_prsd = MatB;
N_3_mid = 3*(Nx*Ny + Nx*Nz + Nx + 2);

MatA_prsd(N_3_mid,:) = [];
MatA_prsd(:,N_3_mid) = [];
MatB_prsd(N_3_mid) = [];
clear("MatA","MatB")
%%

x_mat_raw = mldivide(MatA_prsd,MatB_prsd);
%%

x_mat = zeros(3*N_total_int,1);
for n= 1:3*N_total_int
    if n < N_3_mid
        x_mat(n) = x_mat_raw(n);
    elseif n == N_3_mid
        x_mat(n) = 0.4;
    elseif n > N_3_mid
        x_mat(n) = x_mat_raw(n-1);
    end
end


% x_mat = mldivide(MatA,MatB);

R_mat = zeros(N_total_int,3);

for n = 1:3*N_total_int
    remainder = rem(n,3);
    if (remainder == 0)
        R_mat(n/3 ,3) = x_mat(n);
    elseif (remainder == 1)
        R_mat((n+2)/3,1) = x_mat(n);
    elseif (remainder == 2)
        R_mat((n+1)/3,2) = x_mat(n);
    end

end
