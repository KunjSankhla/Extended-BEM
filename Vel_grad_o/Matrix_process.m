MatA_prsd = MatA;
MatB_prsd = MatB;

i = 0;
for n = 3*N_total+1 : 3*N_total_int
    if rem(n,3) == 0
    else
        MatA_prsd(n-i,:) = [];
        MatA_prsd(:,n-i) = [];
        MatB_prsd(n-i) = [];
        i = i + 1;
    end
end


MatA_prsd(3*N_total+Nx*Ny ,:) = 0;
MatA_prsd(3*N_total+Nx*Ny,3*N_total+1:3*N_total + Nx*Ny) = 1;
MatB_prsd(3*N_total+Nx*Ny) = (Lx*Ly)/(x_diff*y_diff);


%
MatA_prsd(3*N_total,:) = [];
MatA_prsd(:,3*N_total) = [];
MatB_prsd(3*N_total) = [];
%
